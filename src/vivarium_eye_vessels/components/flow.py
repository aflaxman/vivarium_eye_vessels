"""Poiseuille flow on the frozen vessel graph, and shear-driven remodeling.

Hagen-Poiseuille flow through a segment is Ohm's law with conductance
r**4 / L (constants absorbed into the arbitrary pressure units), so the
frozen particle table is a resistor network: parent-child segments and
anastomosis bridges are resistors, artery and vein roots are fixed-pressure
terminals, and Kirchhoff's current law at every junction gives a sparse
linear system for the node pressures. Wall shear stress (what endothelium
actually senses) is then |Q| / r**3 per segment, and remodeling follows
Pries & Secomb: prune terminal segments whose shear falls below a fraction
of the network median, and optionally drift every caliber toward the
median-shear radius. See docs/poiseuille_flow.md for the full story.
"""

from typing import List

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from vivarium import Component
from vivarium.framework.engine import Builder
from vivarium.framework.event import Event

from vivarium_eye_vessels.components.particles import (
    VESSEL_TYPE_ARTERY,
    VESSEL_TYPE_NONE,
    Particle3D,
)

EDGE_COLUMNS = [
    "node_a",
    "node_b",
    "radius",
    "length",
    "conductance",
    "vessel_type",
    "anastomosis",
]


def vessel_edges(pop: pd.DataFrame, min_length: float = 1e-6) -> pd.DataFrame:
    """The frozen vessel graph as edges with Poiseuille conductances.

    One edge per parent-child segment between frozen, calibered particles,
    plus one per anastomosis bridge; conductance is radius**4 / length
    (viscosity constants absorbed into the arbitrary pressure units).
    """
    vessels = pop[pop.frozen & (pop.radius > 0)]
    ids = vessels.index

    def edge_frame(children: pd.DataFrame, others: pd.Series, anastomosis: bool):
        offsets = children[["x", "y", "z"]].to_numpy(float) - pop.loc[
            others, ["x", "y", "z"]
        ].to_numpy(float)
        lengths = np.maximum(np.linalg.norm(offsets, axis=1), min_length)
        radii = children.radius.to_numpy(float)
        return pd.DataFrame(
            {
                "node_a": children.index,
                "node_b": others.to_numpy(),
                "radius": radii,
                "length": lengths,
                "conductance": radii**4 / lengths,
                "vessel_type": children.vessel_type.to_numpy(),
                "anastomosis": anastomosis,
            }
        )

    frames = []
    children = vessels[vessels.parent_id.isin(ids)]
    if not children.empty:
        frames.append(edge_frame(children, children.parent_id, anastomosis=False))
    if "anastomosis_id" in pop.columns:
        joined = vessels[vessels.anastomosis_id.isin(ids)]
        if not joined.empty:
            frames.append(edge_frame(joined, joined.anastomosis_id, anastomosis=True))

    if not frames:
        return pd.DataFrame(columns=EDGE_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def solve_pressures(
    edges: pd.DataFrame,
    boundary: pd.Series,
    leak_conductance: float,
    tissue_pressure: float = 0.0,
) -> pd.Series:
    """Node pressures from Kirchhoff's current law on the vessel graph.

    ``boundary`` maps root nodes to fixed pressures. Every node also leaks
    to a virtual tissue node at ``tissue_pressure`` through
    ``leak_conductance`` — the capillary drainage that happens everywhere
    along real microvessels — which keeps dead-end branches carrying a
    trickle and the system nonsingular even for fragments with no root.
    """
    nodes = pd.Index(sorted(set(edges.node_a) | set(edges.node_b) | set(boundary.index)))
    n = len(nodes)
    positions = pd.Series(np.arange(n), index=nodes)
    i = positions[edges.node_a].to_numpy()
    j = positions[edges.node_b].to_numpy()
    g = edges.conductance.to_numpy(float)

    # Weighted graph Laplacian plus the leak on every diagonal
    rows = np.concatenate([i, j, i, j, np.arange(n)])
    cols = np.concatenate([i, j, j, i, np.arange(n)])
    vals = np.concatenate([g, g, -g, -g, np.full(n, leak_conductance)])
    laplacian = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    rhs = np.full(n, leak_conductance * tissue_pressure)

    fixed_mask = nodes.isin(boundary.index)
    fixed = np.nonzero(fixed_mask)[0]
    free = np.nonzero(~fixed_mask)[0]
    pressures = np.zeros(n)
    pressures[fixed] = boundary.loc[nodes[fixed_mask]].to_numpy(float)
    if len(free):
        # Conductances span many orders of magnitude (r**4), so products can
        # denormalize; flushing them to zero is the correct physics
        with np.errstate(under="ignore"):
            reduced_rhs = rhs[free] - laplacian[free][:, fixed] @ pressures[fixed]
            pressures[free] = spsolve(laplacian[free][:, free].tocsc(), reduced_rhs)
    return pd.Series(pressures, index=nodes)


def edge_flows(edges: pd.DataFrame, pressures: pd.Series) -> pd.DataFrame:
    """Per-edge flow (signed, node_b -> node_a) and wall shear |Q| / r**3."""
    flows = edges.copy()
    pressure_drop = (
        pressures.loc[edges.node_b].to_numpy() - pressures.loc[edges.node_a].to_numpy()
    )
    # Tiny conductance x tiny pressure drop can denormalize; zero is correct
    with np.errstate(under="ignore"):
        flows["flow"] = flows.conductance.to_numpy() * pressure_drop
        flows["shear"] = np.abs(flows.flow) / flows.radius.to_numpy() ** 3
    return flows


class FlowRemodeler(Component):
    """Shear-driven remodeling of the frozen network (roadmap idea 5).

    Every ``remodel_interval`` steps, solve Poiseuille flow on the frozen
    graph (artery roots at ``artery_pressure``, vein roots at
    ``vein_pressure``, a per-node leak to tissue), then:

    - **prune**: terminal segments (degree-1 ends of the graph, so the
      network is never cut in the middle) with shear below
      ``shear_threshold_fraction`` x their own tree's median are recycled —
      the particle's ``unfreeze_time`` is stamped and it returns to the
      free pool as an ordinary wandering particle;
    - **adapt**: every segment's caliber drifts by ``adaptation_rate``
      toward the radius that would put its shear at its own tree's median
      (high-shear trunks thicken, low-shear twigs thin), clipped to
      [``min_radius``, ``max_radius``]. Targets are per tree because
      arteries genuinely run at higher wall shear than veins; a global
      target would erase the artery/vein caliber asymmetry.
    """

    CONFIGURATION_DEFAULTS = {
        "flow_remodeler": {
            "remodel_interval": 25,
            "start_time": "2000-01-01",  # remodel only once the network exists
            "artery_pressure": 1.0,
            "vein_pressure": -1.0,
            "tissue_pressure": 0.0,
            # Per-node leak as a fraction of the median segment conductance
            "leak_fraction": 0.01,
            "shear_threshold_fraction": 0.3,
            "adaptation_rate": 0.0,
            "min_radius": 0.001,
            "max_radius": 0.02,
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return [
            "x",
            "y",
            "z",
            "frozen",
            "freeze_time",
            "unfreeze_time",
            "depth",
            "parent_id",
            "path_id",
            "radius",
            "vessel_type",
            "anastomosis_id",
        ]

    def setup(self, builder: Builder) -> None:
        self.config = builder.configuration.flow_remodeler
        self.clock = builder.time.clock()
        self.start_time = pd.Timestamp(self.config.start_time)
        self.step_count = 0
        self.total_pruned = 0
        self.particles = builder.components.get_components_by_type(Particle3D)[0]

    def on_time_step(self, event: Event) -> None:
        self.step_count += 1
        if self.step_count % self.config.remodel_interval:
            return
        if self.clock() < self.start_time:
            return
        pop = self.population_view.get(event.index, self.required_attributes)
        flows = self.solve_network(pop)
        if flows is None:
            return
        self.remodel(pop, flows)

    def solve_network(self, pop: pd.DataFrame) -> pd.DataFrame | None:
        """Poiseuille flow and shear per frozen segment, or None if unsolvable."""
        edges = vessel_edges(pop)
        if edges.empty:
            return None
        roots = pop[
            pop.frozen & (pop.parent_id < 0) & (pop.path_id >= 0) & (pop.vessel_type > 0)
        ]
        if roots.empty:
            return None
        boundary = pd.Series(
            np.where(
                roots.vessel_type == VESSEL_TYPE_ARTERY,
                float(self.config.artery_pressure),
                float(self.config.vein_pressure),
            ),
            index=roots.index,
        )
        leak = float(self.config.leak_fraction) * float(edges.conductance.median())
        pressures = solve_pressures(edges, boundary, leak, float(self.config.tissue_pressure))
        return edge_flows(edges, pressures)

    def graph_degrees(self, pop: pd.DataFrame) -> pd.Series:
        """Connection count per particle over every current or future edge.

        Counts the parent link, links from children (frozen or still
        active), the particle's own anastomosis bridge, and bridges from
        other particles onto it — so a "terminal" segment (degree 1) is
        truly the end of a branch and pruning never cuts the graph.
        """
        on_graph = pop[(pop.path_id >= 0) | pop.frozen]
        degrees = pd.Series(0, index=pop.index)
        has_parent = on_graph[on_graph.parent_id.isin(pop.index)]
        degrees[has_parent.index] += 1
        degrees = degrees.add(has_parent.parent_id.value_counts(), fill_value=0)
        bridged = on_graph[on_graph.anastomosis_id.isin(pop.index)]
        degrees[bridged.index] += 1
        degrees = degrees.add(bridged.anastomosis_id.value_counts(), fill_value=0)
        return degrees

    def remodel(self, pop: pd.DataFrame, flows: pd.DataFrame) -> None:
        positive = flows[flows.shear > 0]
        if positive.empty:
            return
        # Each tree remodels toward its own median shear: arteries genuinely
        # run at higher wall shear than veins, so a global target would erase
        # the artery/vein caliber asymmetry (the clinical AVR)
        type_medians = positive.groupby("vessel_type").shear.median()
        median_shear = flows.vessel_type.map(type_medians).fillna(
            float(positive.shear.median())
        )

        # --- Prune low-shear terminal segments ---
        thresholds = float(self.config.shear_threshold_fraction) * median_shear
        degrees = self.graph_degrees(pop)
        low_shear_children = flows.node_a[flows.shear < thresholds]
        candidates = pop.loc[low_shear_children]
        candidates = candidates[
            candidates.frozen & (candidates.parent_id >= 0) & (degrees[candidates.index] == 1)
        ]
        pruned = candidates.index
        if not pruned.empty:
            self.total_pruned += len(pruned)
            recycled = pd.DataFrame(
                {
                    "frozen": False,
                    "freeze_time": pd.NaT,
                    "unfreeze_time": self.clock(),
                    "path_id": -1,
                    "parent_id": -1,
                    "depth": -1,
                    "radius": 0.0,
                    "vessel_type": VESSEL_TYPE_NONE,
                    "anastomosis_id": -1,
                },
                index=pruned,
            )
            self.particles.update_particles(recycled)

        # --- Adapt calibers toward the median-shear radius ---
        rate = float(self.config.adaptation_rate)
        if rate <= 0:
            return
        # The depth-0 arcades are boundary conditions — their calibers are set
        # upstream by the central retinal artery and vein (and carry the
        # seeded A:V ratio, the clinical biomarker) — so they don't adapt
        arcade = pop.loc[flows.node_a, "depth"].to_numpy() == 0
        keep = ~flows.anastomosis & ~flows.node_a.isin(pruned) & ~arcade
        segments = flows[keep]
        if segments.empty:
            return
        with np.errstate(divide="ignore", under="ignore"):
            factors = np.power(
                np.clip(
                    segments.shear.to_numpy() / median_shear[keep].to_numpy(), 1e-6, None
                ),
                rate,
            )
        new_radii = np.clip(
            segments.radius.to_numpy() * factors,
            float(self.config.min_radius),
            float(self.config.max_radius),
        )
        self.particles.update_particles(
            pd.DataFrame({"radius": new_radii}, index=pd.Index(segments.node_a))
        )

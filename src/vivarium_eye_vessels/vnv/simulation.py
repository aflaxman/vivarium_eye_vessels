"""Headless construction and execution of model specifications.

The model specifications in this repository include the pygame visualizer,
which is unhelpful for batch V&V runs. This module builds an
:class:`~vivarium.InteractiveContext` from a model specification YAML while
skipping visualization components.
"""

import copy
import importlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml
from vivarium import InteractiveContext

SKIPPED_MODULES = ("visualizer",)

TREE_ATTRIBUTES = [
    "x",
    "y",
    "z",
    "vx",
    "vy",
    "vz",
    "frozen",
    "path_id",
    "parent_id",
    "depth",
    "radius",
    "vessel_type",
    "anastomosis_id",
    "layer_id",
]


def build_headless_simulation(spec_path: str | Path) -> InteractiveContext:
    """Build an InteractiveContext from a model spec, skipping the visualizer."""
    with open(spec_path) as f:
        spec = yaml.safe_load(f)

    components = []
    for package, entries in spec.get("components", {}).items():
        for entry in entries:
            entry = entry.strip()
            if not entry.endswith("()"):
                raise NotImplementedError(
                    f"Only no-argument component constructors are supported, got: {entry}"
                )
            module_name, class_name = entry[:-2].rsplit(".", 1)
            if module_name.split(".")[0] in SKIPPED_MODULES:
                continue
            module = importlib.import_module(f"{package}.{module_name}")
            components.append(getattr(module, class_name)())

    configuration = spec.get("configuration", {})
    return InteractiveContext(components=components, configuration=configuration)


def with_seed(spec: dict, seed: int) -> dict:
    """A deep copy of ``spec`` with its random seed replaced."""
    candidate = copy.deepcopy(spec)
    candidate["configuration"]["randomness"]["random_seed"] = seed
    return candidate


def build_from_spec(spec: dict, spec_path: str | Path) -> InteractiveContext:
    """Write ``spec`` to ``spec_path`` (a record of the run) and build it headless."""
    with open(spec_path, "w") as f:
        yaml.safe_dump(spec, f)
    return build_headless_simulation(spec_path)


def get_ellipsoid_bounds(sim: InteractiveContext) -> tuple[float, float]:
    """Get the (a, b) semi-axes of the containment ellipsoid, if configured."""
    a, b, _ = get_ellipsoid_semi_axes(sim)
    return a, b


def get_ellipsoid_semi_axes(sim: InteractiveContext) -> tuple[float, float, float]:
    """Get the (a, b, c) semi-axes of the containment ellipsoid, if configured."""
    try:
        config = sim.configuration.ellipsoid_containment
        return float(config.a), float(config.b), float(config.c)
    except AttributeError:
        return 2.0, 2.0, 2.0


def get_perfusion_params(sim: InteractiveContext) -> tuple[float, float]:
    """The (site_spacing, perfusion_radius) the PerfusionDemand component runs with.

    The perfused-fraction metric mirrors that component's site lattice, so
    it must read the same numbers the model does rather than assume them.
    """
    try:
        config = sim.configuration.perfusion_demand
        return float(config.site_spacing), float(config.perfusion_radius)
    except AttributeError:
        return 0.1, 0.15


def get_disc_center(sim: InteractiveContext) -> tuple[float, float]:
    """The (x, y) of the optic disc: the center of the root circle."""
    try:
        center = sim.configuration.particles.initial_circle.center
        return float(center[0]), float(center[1])
    except (AttributeError, IndexError, TypeError):
        return 0.0, 0.0


@dataclass(frozen=True)
class Geometry:
    """The spatial facts of a model the V&V metrics have to agree with.

    Read from the built simulation (:func:`get_geometry`) rather than
    assumed, so a spec change moves the metrics with it.
    """

    semi_axes: tuple[float, float, float]  # containment ellipsoid (a, b, c)
    perfusion: tuple[float, float]  # PerfusionDemand (site_spacing, perfusion_radius)
    disc_center: tuple[float, float]  # optic disc (x, y): the root circle's center

    @property
    def bounds(self) -> tuple[float, float]:
        """The (a, b) semi-axes the x-y raster spans."""
        return self.semi_axes[0], self.semi_axes[1]


def get_geometry(sim: InteractiveContext) -> Geometry:
    return Geometry(
        semi_axes=get_ellipsoid_semi_axes(sim),
        perfusion=get_perfusion_params(sim),
        disc_center=get_disc_center(sim),
    )


def get_network(sim: InteractiveContext) -> pd.DataFrame:
    """Get the particle table attributes needed for network analysis."""
    return sim.get_population(TREE_ATTRIBUTES)


def tree_edges(pop: pd.DataFrame) -> pd.DataFrame:
    """Get vessel segments as a frame of parent/child coordinate pairs.

    Segments connect each particle with a valid parent to that parent's
    position; particles whose parent has left the table are skipped.
    Anastomoses (tips fused onto the other tree) contribute one extra
    segment each, from the tip to its ``anastomosis_id`` target, flagged
    in the ``anastomosis`` column.
    """

    def edge_frame(children: pd.DataFrame, other_ids: pd.Series, anastomosis: bool):
        others = pop.loc[other_ids]
        return pd.DataFrame(
            {
                "x0": others.x.values,
                "y0": others.y.values,
                "z0": others.z.values,
                "x1": children.x.values,
                "y1": children.y.values,
                "z1": children.z.values,
                "child": children.index.values,
                "parent": other_ids.values,
                "radius": children.radius.values,
                "vessel_type": children.vessel_type.values,
                "layer_id": children.layer_id.values,
                "anastomosis": anastomosis,
            }
        )

    children = pop[pop.parent_id >= 0]
    children = children[children.parent_id.isin(pop.index)]
    frames = [edge_frame(children, children.parent_id, anastomosis=False)]

    if "anastomosis_id" in pop.columns:
        joined = pop[pop.anastomosis_id >= 0]
        joined = joined[joined.anastomosis_id.isin(pop.index)]
        if not joined.empty:
            frames.append(edge_frame(joined, joined.anastomosis_id, anastomosis=True))

    return pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]


def run_steps(sim: InteractiveContext, n_steps: int, callback=None, every: int = 1) -> None:
    """Step the simulation, optionally invoking a callback every ``every`` steps."""
    for step in range(n_steps):
        sim.step()
        if callback is not None and (step + 1) % every == 0:
            callback(step + 1)

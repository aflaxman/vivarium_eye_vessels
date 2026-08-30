"""Tests for Poiseuille flow and shear-driven remodeling (roadmap idea 5)."""

import copy

import numpy as np
import pandas as pd
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import (
    EllipsoidContainment,
    FrozenRepulsion,
    PerfusionDemand,
)
from vivarium_eye_vessels.components.flow import (
    FlowRemodeler,
    edge_flows,
    solve_pressures,
    vessel_edges,
)
from vivarium_eye_vessels.components.particles import (
    VESSEL_TYPE_ARTERY,
    VESSEL_TYPE_NONE,
    VESSEL_TYPE_VEIN,
    Particle3D,
    PathAnastomosis,
    PathExtinction,
    PathFreezer,
    PathSplitter,
)
from vivarium_eye_vessels.vnv import simulation


def straight_pipe_population() -> pd.DataFrame:
    """Artery root -> middle -> vein-joined tip: one artery-to-vein pipe."""
    return pd.DataFrame(
        {
            "x": [0.0, 1.0, 2.0, 3.0],
            "y": [0.0, 0.0, 0.0, 0.0],
            "z": [0.0, 0.0, 0.0, 0.0],
            "frozen": [True, True, True, True],
            "parent_id": [-1, 0, 1, -1],
            "path_id": [0, 0, 0, 1],
            "depth": [0, 0, 0, 0],
            "radius": [0.01, 0.01, 0.01, 0.01],
            "vessel_type": [
                VESSEL_TYPE_ARTERY,
                VESSEL_TYPE_ARTERY,
                VESSEL_TYPE_ARTERY,
                VESSEL_TYPE_VEIN,
            ],
            "anastomosis_id": [-1, -1, 3, -1],
        },
        index=[0, 1, 2, 3],
    )


def test_single_pipe_flow_matches_hand_calculation():
    pop = straight_pipe_population()
    edges = vessel_edges(pop)
    assert len(edges) == 3  # two segments plus one anastomosis bridge
    boundary = pd.Series({0: 1.0, 3: -1.0})
    # Negligible leak: pure series circuit, conductance g = r^4 / L each
    pressures = solve_pressures(edges, boundary, leak_conductance=1e-15)
    flows = edge_flows(edges, pressures)

    g = 0.01**4 / 1.0
    expected_flow = (1.0 - (-1.0)) / 3 * g  # three equal resistors in series
    np.testing.assert_allclose(np.abs(flows.flow), expected_flow, rtol=1e-6)
    np.testing.assert_allclose(flows.shear, expected_flow / 0.01**3, rtol=1e-6)
    # Pressure falls monotonically from artery to vein
    assert pressures[0] > pressures[1] > pressures[2] > pressures[3]


def test_parallel_pipes_split_flow_by_fourth_power_of_radius():
    pop = pd.DataFrame(
        {
            "x": [0.0, 1.0, 1.0],
            "y": [0.0, 1.0, -1.0],
            "z": [0.0, 0.0, 0.0],
            "frozen": [True, True, True],
            "parent_id": [-1, 0, 0],
            "path_id": [0, 1, 2],
            "depth": [0, 1, 1],
            "radius": [0.02, 0.01, 0.02],
            "vessel_type": [VESSEL_TYPE_ARTERY] * 3,
            "anastomosis_id": [-1, -1, -1],
        },
        index=[0, 1, 2],
    )
    edges = vessel_edges(pop)
    boundary = pd.Series({0: 1.0, 1: 0.0, 2: 0.0})
    pressures = solve_pressures(edges, boundary, leak_conductance=1e-15)
    flows = edge_flows(edges, pressures).set_index("node_a")

    # Same length and pressure drop, so flow ratio is the conductance ratio 2^4
    np.testing.assert_allclose(flows.loc[2].flow / flows.loc[1].flow, 16.0, rtol=1e-6)


def test_leak_gives_dead_ends_nonzero_flow():
    pop = straight_pipe_population().drop(index=3)
    pop["anastomosis_id"] = -1  # no bridge: pure dead-end artery branch
    edges = vessel_edges(pop)
    boundary = pd.Series({0: 1.0})
    pressures = solve_pressures(edges, boundary, leak_conductance=1e-9, tissue_pressure=0.0)
    flows = edge_flows(edges, pressures)
    assert (np.abs(flows.flow) > 0).all()
    # Flow accumulates toward the root: the root-side segment drains more leaks
    by_child = flows.set_index("node_a")
    assert np.abs(by_child.loc[1].flow) > np.abs(by_child.loc[2].flow)


def make_remodeling_simulation(**overrides) -> tuple[InteractiveContext, FlowRemodeler]:
    config = {
        "randomness": {"random_seed": 42},
        "time": {
            "start": {"year": 2025, "month": 1, "day": 1},
            "end": {"year": 2025, "month": 6, "day": 1},
            "step_size": 0.05,
        },
        "population": {"population_size": 300},
        "particles": {
            "overall_max_velocity_change": 0.2,
            "initial_velocity_range": [-0.2, 0.2],
            "terminal_velocity": 0.15,
            "initial_circle": {"center": [1.0, 0.0, 0.0], "radius": 0.05, "n_vessels": 4},
            "root_radius": 0.02,
            "artery_caliber_ratio": 0.67,
        },
        "path_freezer": {"freeze_interval": 3, "radius_taper": 0.99},
        "path_splitter": {
            "split_interval": 15,
            "split_angle": 60,
            "split_probability": 0.9,
            "max_depth": 4,
            "murray_exponent": 3.0,
            "flow_asymmetry": 0.15,
            "min_radius": 0.002,
            "caliber_cadence_exponent": 0.0,
        },
        "path_extinction": {"force_threshold": 1.2},
        "path_anastomosis": {
            "capture_radius": 0.05,
            "max_tip_radius": 0.01,
            "max_target_radius": 0.01,
            "probability": 0.9,
        },
        "flow_remodeler": {
            "remodel_interval": 20,
            "start_time": "2025-01-03",
            "artery_pressure": 1.0,
            "vein_pressure": -1.0,
            "tissue_pressure": 0.0,
            "leak_fraction": 0.01,
            "shear_threshold_fraction": 0.5,
            "adaptation_rate": 0.05,
            "min_radius": 0.001,
            "max_radius": 0.02,
        },
        "frozen_repulsion": {
            "force_type": "hookean",
            "spring_constant": 1.5,
            "interaction_radius": 0.15,
            "delay": 1,
            "cross_type_factor": 0.25,
        },
        "ellipsoid_containment": {
            "a": 2,
            "b": 2,
            "c": 0.2,
            "force_type": "hookean",
            "spring_constant": 3,
        },
        "perfusion_demand": {
            "site_spacing": 0.1,
            "perfusion_radius": 0.15,
            "influence_radius": 2.0,
            "magnitude": 0.3,
        },
    }
    for section, values in overrides.items():
        config[section].update(values)
    remodeler = FlowRemodeler()
    sim = InteractiveContext(
        components=[
            Particle3D(),
            PathFreezer(),
            PathExtinction(),
            PathSplitter(),
            PathAnastomosis(),
            remodeler,
            EllipsoidContainment(),
            FrozenRepulsion(),
            PerfusionDemand(),
        ],
        configuration=copy.deepcopy(config),
    )
    return sim, remodeler


def test_remodeler_prunes_and_recycles_without_cutting_the_graph():
    sim, remodeler = make_remodeling_simulation()

    # Step to just past the first pruning pass. The remodeler runs after the
    # growth components within a step, so freshly recycled particles cannot
    # have been re-used yet when we look.
    for _ in range(250):
        pruned_before = remodeler.total_pruned
        sim.step()
        if remodeler.total_pruned > pruned_before:
            break
    assert remodeler.total_pruned > 0, "no segments were pruned"

    pop = sim.get_population(
        [
            "frozen",
            "path_id",
            "parent_id",
            "radius",
            "vessel_type",
            "anastomosis_id",
            "unfreeze_time",
        ]
    )

    # Recycled particles carry the unfreeze_time stamp and are free again
    recycled = pop[pop.unfreeze_time.notna() & (pop.path_id < 0) & ~pop.frozen]
    assert len(recycled) > 0
    assert (recycled.vessel_type == VESSEL_TYPE_NONE).all()
    assert (recycled.radius == 0).all()

    # Graph integrity: no vessel particle points at a recycled particle
    vessels = pop[(pop.path_id >= 0) | pop.frozen]
    vessels = vessels[vessels.radius > 0]
    parents_ok = ~vessels.parent_id.isin(recycled.index)
    bridges_ok = ~vessels.anastomosis_id.isin(recycled.index)
    assert parents_ok.all(), "a pruned particle is still someone's parent"
    assert bridges_ok.all(), "a pruned particle is still an anastomosis target"


def test_adaptation_widens_shear_spread_of_calibers():
    """With adaptation on, trunk and twig calibers separate further.

    Deterministic seeds make this an exact comparison, not a statistical one.
    """

    def radius_spread(adaptation_rate: float) -> float:
        sim, _ = make_remodeling_simulation(
            flow_remodeler={"adaptation_rate": adaptation_rate}
        )
        simulation.run_steps(sim, 200)
        pop = sim.get_population(["frozen", "radius", "path_id"])
        radii = pop[pop.frozen & (pop.radius > 0)].radius
        return float(radii.quantile(0.9) / radii.quantile(0.1))

    assert radius_spread(0.2) > radius_spread(0.0)

"""Tests for the growth-reliability mechanisms: balanced arterial inflow
(no rich-get-richer arcade) and the crowding gate on splits."""

import copy

import numpy as np
import pandas as pd
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import EllipsoidContainment
from vivarium_eye_vessels.components.flow import (
    edge_flows,
    root_outflows,
    solve_pressures,
    vessel_edges,
)
from vivarium_eye_vessels.components.particles import (
    VESSEL_TYPE_ARTERY,
    VESSEL_TYPE_VEIN,
    Particle3D,
    PathFreezer,
    PathSplitter,
)
from vivarium_eye_vessels.vnv import simulation


def two_arcade_population() -> pd.DataFrame:
    """Two artery arcades of very different caliber, draining into one vein."""
    return pd.DataFrame(
        {
            "x": [0.0, 0.0, 1.0, 1.0, 2.0, 2.0],
            "y": [0.0, 1.0, 0.0, 1.0, 0.5, 0.5],
            "z": [0.0] * 6,
            "frozen": [True] * 6,
            "parent_id": [-1, -1, 0, 1, -1, 4],
            "path_id": [0, 1, 0, 1, 2, 2],
            "depth": [0] * 6,
            # Arcade A (nodes 0->2) is twice the caliber of arcade B (1->3),
            # a 16x conductance advantage
            "radius": [0.02, 0.01, 0.02, 0.01, 0.02, 0.02],
            "vessel_type": [
                VESSEL_TYPE_ARTERY,
                VESSEL_TYPE_ARTERY,
                VESSEL_TYPE_ARTERY,
                VESSEL_TYPE_ARTERY,
                VESSEL_TYPE_VEIN,
                VESSEL_TYPE_VEIN,
            ],
            # Both arcade tips fuse onto the vein branch
            "anastomosis_id": [-1, -1, 5, 5, -1, -1],
        },
        index=range(6),
    )


def test_dirichlet_roots_let_the_wide_arcade_dominate():
    edges = vessel_edges(two_arcade_population())
    boundary = pd.Series({0: 1.0, 1: 1.0, 4: -1.0})
    pressures = solve_pressures(edges, boundary, leak_conductance=1e-15)
    outflows = root_outflows(edge_flows(edges, pressures), pd.Index([0, 1]))
    assert outflows[0] > 5 * outflows[1], "fixed pressures should favor the wide arcade"


def test_balanced_injections_split_inflow_equally():
    edges = vessel_edges(two_arcade_population())
    injections = pd.Series({0: 0.5, 1: 0.5})
    pressures = solve_pressures(
        edges, pd.Series({4: -1.0}), leak_conductance=1e-15, injections=injections
    )
    outflows = root_outflows(edge_flows(edges, pressures), pd.Index([0, 1]))
    np.testing.assert_allclose(outflows[0], 0.5, rtol=1e-6)
    np.testing.assert_allclose(outflows[1], 0.5, rtol=1e-6)
    # The narrow arcade now runs at much higher shear than the wide one, so
    # shear-driven adaptation thickens the under-built arcade
    flows = edge_flows(edges, pressures)
    tree = flows[~flows.anastomosis].set_index("node_a")
    assert tree.loc[3].shear > tree.loc[2].shear


CONFIGURATION = {
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
    },
    "path_freezer": {"freeze_interval": 3, "radius_taper": 0.999},
    "path_splitter": {
        "split_interval": 20,
        "split_angle": 60,
        "split_probability": 0.9,
        "max_depth": 4,
        "murray_exponent": 3.0,
        "flow_asymmetry": 0.15,
        "min_radius": 0.002,
        "crowding_radius": 0.06,
        "max_crowding": 1,
    },
    "ellipsoid_containment": {
        "a": 2,
        "b": 2,
        "c": 0.2,
        "force_type": "hookean",
        "spring_constant": 3,
    },
}


def test_crowding_gate_drops_saturated_tips_and_keeps_the_front():
    splitter = PathSplitter()
    sim = InteractiveContext(
        components=[Particle3D(), PathFreezer(), splitter, EllipsoidContainment()],
        configuration=copy.deepcopy(CONFIGURATION),
    )
    simulation.run_steps(sim, 30)
    frozen = splitter.freezer.frozen_particles()
    assert frozen is not None and len(frozen) > 0

    # One candidate sits exactly on a frozen particle (saturated at limit 1),
    # the other far outside the network (empty surroundings)
    crowded_position = frozen[["x", "y", "z"]].iloc[0]
    candidates = pd.DataFrame(
        {
            "x": [crowded_position.x, 10.0],
            "y": [crowded_position.y, 10.0],
            "z": [crowded_position.z, 10.0],
        },
        index=[900, 901],
    )
    kept = splitter.uncrowded(candidates, pd.Index([900, 901]))
    assert list(kept) == [901]

    # Limit 0 disables the gate entirely
    config = copy.deepcopy(CONFIGURATION)
    config["path_splitter"]["max_crowding"] = 0
    ungated = PathSplitter()
    InteractiveContext(
        components=[Particle3D(), PathFreezer(), ungated, EllipsoidContainment()],
        configuration=config,
    )
    assert list(ungated.uncrowded(candidates, pd.Index([900, 901]))) == [900, 901]

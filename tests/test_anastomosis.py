"""Tests for anastomosis: capillary loops between the trees (roadmap idea 4)."""

import copy

import numpy as np
import pandas as pd
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import (
    EllipsoidContainment,
    FrozenRepulsion,
    PerfusionDemand,
)
from vivarium_eye_vessels.components.particles import (
    VESSEL_TYPE_ARTERY,
    VESSEL_TYPE_VEIN,
    Particle3D,
    PathAnastomosis,
    PathExtinction,
    PathFreezer,
    PathSplitter,
    anastomosis_targets,
)
from vivarium_eye_vessels.vnv import metrics, simulation

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


def test_anastomosis_targets_picks_nearest_opposite_capillary():
    tips = pd.DataFrame(
        {
            "x": [0.0],
            "y": [0.0],
            "z": [0.0],
            "vessel_type": [VESSEL_TYPE_ARTERY],
            "radius": [0.003],
        },
        index=[10],
    )
    frozen = pd.DataFrame(
        {
            "x": [0.02, 0.01, 0.015, 0.01],
            "y": [0.0, 0.0, 0.0, 0.0],
            "z": [0.0, 0.0, 0.0, 0.0],
            # same tree, wide trunk, nearest vein capillary, farther vein
            "vessel_type": [
                VESSEL_TYPE_ARTERY,
                VESSEL_TYPE_VEIN,
                VESSEL_TYPE_VEIN,
                VESSEL_TYPE_VEIN,
            ],
            "radius": [0.003, 0.02, 0.003, 0.003],
        },
        index=[100, 101, 102, 103],
    )
    targets = anastomosis_targets(tips, frozen, [[0, 1, 2, 3]], max_target_radius=0.01)
    assert list(targets.index) == [10]
    assert targets[10] == 103  # nearest vein with capillary caliber; trunk 101 excluded


def test_anastomosis_targets_empty_without_opposite_type():
    tips = pd.DataFrame(
        {
            "x": [0.0],
            "y": [0.0],
            "z": [0.0],
            "vessel_type": [VESSEL_TYPE_ARTERY],
            "radius": [0.003],
        },
        index=[10],
    )
    frozen = pd.DataFrame(
        {
            "x": [0.01],
            "y": [0.0],
            "z": [0.0],
            "vessel_type": [VESSEL_TYPE_ARTERY],
            "radius": [0.003],
        },
        index=[100],
    )
    targets = anastomosis_targets(tips, frozen, [[0]], max_target_radius=0.01)
    assert targets.empty


def test_graph_cycles_counts_anastomosis_loops():
    # Two chains bridged twice form one loop; bridged once they are still a tree
    pop = pd.DataFrame(
        {
            "parent_id": [-1, 0, -1, 2],
            "anastomosis_id": [-1, 2, -1, 0],
        },
        index=[0, 1, 2, 3],
    )
    assert metrics.graph_cycles(pop) == 1

    # A single bridge merges the trees without creating a cycle
    pop["anastomosis_id"] = [-1, -1, -1, 0]
    assert metrics.graph_cycles(pop) == 0

    # No bridges: a forest has zero cycles
    pop["anastomosis_id"] = -1
    assert metrics.graph_cycles(pop) == 0


def test_simulation_forms_anastomoses_and_cycles():
    sim = InteractiveContext(
        components=[
            Particle3D(),
            PathFreezer(),
            PathExtinction(),
            PathSplitter(),
            PathAnastomosis(),
            EllipsoidContainment(),
            FrozenRepulsion(),
            PerfusionDemand(),
        ],
        configuration=copy.deepcopy(CONFIGURATION),
    )
    simulation.run_steps(sim, 250)
    pop = sim.get_population(
        ["frozen", "path_id", "radius", "vessel_type", "anastomosis_id", "parent_id"]
    )

    joined = pop[pop.anastomosis_id >= 0]
    assert len(joined) > 0, "no anastomoses formed"

    # Joined tips froze and their paths ended
    assert joined.frozen.all()
    assert (joined.path_id == -1).all()

    # Every join connects opposite trees at capillary calibers
    targets = pop.loc[joined.anastomosis_id]
    assert (targets.vessel_type.values != joined.vessel_type.values).all()
    assert (targets.vessel_type.values > 0).all()
    max_target = CONFIGURATION["path_anastomosis"]["max_target_radius"]
    assert (joined.radius <= CONFIGURATION["path_anastomosis"]["max_tip_radius"]).all()
    assert (targets.radius <= max_target).all()

    # The joins turn the trees into a graph with cycles. A join that merges
    # two previously separate trees adds no cycle, and there are 4 root
    # trees, so at most 3 joins are merges; every other join closes a loop.
    assert metrics.graph_cycles(pop) >= len(joined) - 3

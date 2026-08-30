"""Tests for caliber-dependent branching cadence (roadmap idea 1b)."""

import copy

import numpy as np
import pandas as pd
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import EllipsoidContainment
from vivarium_eye_vessels.components.particles import (
    Particle3D,
    PathFreezer,
    PathSplitter,
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
        "caliber_cadence_exponent": 0.0,
    },
    "ellipsoid_containment": {
        "a": 2,
        "b": 2,
        "c": 0.2,
        "force_type": "hookean",
        "spring_constant": 3,
    },
}


def make_simulation(cadence_exponent: float) -> tuple[InteractiveContext, PathSplitter]:
    config = copy.deepcopy(CONFIGURATION)
    config["path_splitter"]["caliber_cadence_exponent"] = cadence_exponent
    splitter = PathSplitter()
    sim = InteractiveContext(
        components=[Particle3D(), PathFreezer(), splitter, EllipsoidContainment()],
        configuration=config,
    )
    return sim, splitter


def test_zero_exponent_gives_uniform_probability():
    _, splitter = make_simulation(0.0)
    active = pd.DataFrame({"radius": [0.02, 0.005, 0.002, 0.0]}, index=[3, 5, 8, 13])
    probabilities = splitter.split_probabilities(active)
    assert list(probabilities.index) == [3, 5, 8, 13]
    np.testing.assert_allclose(probabilities, 0.9)


def test_probabilities_scale_with_caliber():
    _, splitter = make_simulation(1.0)
    active = pd.DataFrame({"radius": [0.002, 0.004, 0.02, 0.0, 0.001]}, index=[1, 2, 3, 4, 5])
    probabilities = splitter.split_probabilities(active)

    # Base probability applies at the min_radius caliber floor
    np.testing.assert_allclose(probabilities[1], 0.9)
    # Wider tips split less often: (min_radius / radius) ** exponent
    np.testing.assert_allclose(probabilities[2], 0.45)
    np.testing.assert_allclose(probabilities[3], 0.09)
    # Uncalibered tips keep the base probability
    np.testing.assert_allclose(probabilities[4], 0.9)
    # Below the floor the factor exceeds 1 and is clipped to a valid probability
    np.testing.assert_allclose(probabilities[5], 1.0)


def test_cadence_lengthens_trunk_segments():
    """Wide-caliber tips branching less often means longer tree segments.

    Deterministic seeds make this an exact comparison, not a statistical one.
    """

    def median_segment_length(cadence_exponent: float) -> float:
        sim, _ = make_simulation(cadence_exponent)
        simulation.run_steps(sim, 250)
        pop = sim.get_population(["x", "y", "z", "frozen", "path_id", "parent_id"])
        assert pop.frozen.sum() > 100, "network failed to grow"
        lengths = metrics.tree_segment_lengths(pop)
        assert len(lengths) > 0
        return float(np.median(lengths))

    assert median_segment_length(1.5) > median_segment_length(0.0)


def test_tree_segment_lengths_on_known_tree():
    """A hand-built Y tree: trunk of two edges, then two single-edge daughters."""
    pop = pd.DataFrame(
        {
            "x": [0.0, 1.0, 2.0, 2.0, 2.0],
            "y": [0.0, 0.0, 0.0, 1.0, -3.0],
            "z": [0.0, 0.0, 0.0, 0.0, 0.0],
            "path_id": [0, 0, 0, 1, 2],
            "parent_id": [-1, 0, 1, 2, 2],
        },
        index=[0, 1, 2, 3, 4],
    )
    lengths = metrics.tree_segment_lengths(pop)
    np.testing.assert_allclose(sorted(lengths), [1.0, 2.0, 3.0])

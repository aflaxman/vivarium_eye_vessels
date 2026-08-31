"""Tests for Ornstein-Uhlenbeck steering and controllable tortuosity (idea 7)."""

import copy

import numpy as np
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
        "noise_persistence_time": 0.0,
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


def make_simulation(persistence_time: float) -> InteractiveContext:
    config = copy.deepcopy(CONFIGURATION)
    config["particles"]["noise_persistence_time"] = persistence_time
    return InteractiveContext(
        components=[Particle3D(), PathFreezer(), PathSplitter(), EllipsoidContainment()],
        configuration=config,
    )


def test_legacy_mode_leaves_wander_state_untouched():
    sim = make_simulation(0.0)
    simulation.run_steps(sim, 30)
    pop = sim.get_population(["wx", "wy", "wz"])
    assert (pop[["wx", "wy", "wz"]] == 0.0).all().all()


def test_ou_steering_has_configured_statistics():
    """Stationary spread and lag-1 autocorrelation match the AR(1) design."""
    sim = make_simulation(0.5)
    theta = 0.05 / 0.5
    for _ in range(60):
        sim.step()
    # Watch a fixed set of free wanderers: the population grows over time and
    # frozen particles stop updating their steering state
    pop = sim.get_population(["wx", "frozen", "path_id"])
    watched = pop[~pop.frozen & (pop.path_id < 0)].index
    previous = pop.wx.loc[watched].to_numpy()
    correlations = []
    for _ in range(30):
        sim.step()
        current = sim.get_population(["wx"]).wx.loc[watched].to_numpy()
        correlations.append(np.corrcoef(previous, current)[0, 1])
        previous = current

    assert (previous != 0).any()
    # Stationary sd = max_velocity_change / sqrt(3) x scale_x = 0.2/sqrt(3) x 2
    np.testing.assert_allclose(np.std(previous), 0.2 / np.sqrt(3) * 2.0, rtol=0.25)
    np.testing.assert_allclose(np.mean(correlations), 1 - theta, atol=0.08)


def test_persistence_reduces_tortuosity():
    """Longer steering memory means straighter vessels.

    Deterministic seeds make this an exact comparison, not a statistical one.
    """

    def median_tortuosity(persistence_time: float) -> float:
        sim = make_simulation(persistence_time)
        simulation.run_steps(sim, 250)
        pop = sim.get_population(["x", "y", "z", "frozen", "path_id", "parent_id", "depth"])
        assert pop.frozen.sum() > 100, "network failed to grow"
        ratios = metrics.path_tortuosity(pop)
        assert len(ratios) > 0
        return float(np.median(ratios))

    wiggly = median_tortuosity(0.05)  # correlation time of a single step
    smooth = median_tortuosity(2.0)
    assert smooth < wiggly

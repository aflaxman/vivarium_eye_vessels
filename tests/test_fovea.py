"""Tests for the foveal avascular zone: caliber-aware exclusion and demand-free fovea."""

import copy

import numpy as np
import pandas as pd
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import (
    CylinderExclusion,
    EllipsoidContainment,
    PerfusionDemand,
)
from vivarium_eye_vessels.components.particles import Particle3D, PathFreezer

CONFIGURATION = {
    "randomness": {"random_seed": 1},
    "time": {
        "start": {"year": 2025, "month": 1, "day": 1},
        "end": {"year": 2025, "month": 2, "day": 1},
        "step_size": 0.05,
    },
    "population": {"population_size": 50},
    "particles": {
        "initial_circle": {"center": [1.0, 0.0, 0.0], "radius": 0.05, "n_vessels": 4},
        "root_radius": 0.02,
    },
    "ellipsoid_containment": {"a": 2, "b": 2, "c": 0.1, "force_type": "hookean"},
    "perfusion_demand": {"site_spacing": 0.05, "perfusion_radius": 0.15},
    "cylinder_exclusion": {
        "radius": 0.1,
        "wide_radius": 0.2,
        "wide_min_radius": 0.0025,
        "center": [0.0, 0.0, 0.0],
        "direction": [0.0, 0.0, 1.0],
        "force_type": "hookean",
        "spring_constant": 10,
    },
}


def build(**overrides) -> tuple[CylinderExclusion, PerfusionDemand]:
    config = copy.deepcopy(CONFIGURATION)
    config["cylinder_exclusion"].update(overrides)
    exclusion, demand = CylinderExclusion(), PerfusionDemand()
    InteractiveContext(
        components=[Particle3D(), PathFreezer(), EllipsoidContainment(), exclusion, demand],
        configuration=config,
    )
    return exclusion, demand


def tips(x, radii) -> pd.DataFrame:
    return pd.DataFrame({"x": x, "y": 0.0, "z": 0.0, "frozen": False, "radius": radii})


def test_wide_tips_are_held_farther_from_the_fovea_than_capillaries():
    exclusion, _ = build()
    # Two tips 0.15 from the fovea: inside the wide radius, outside the capillary one
    forces = exclusion.calculate_forces_vectorized(tips([0.15, 0.15], [0.002, 0.005]))
    np.testing.assert_allclose(forces[0], 0.0)  # a capillary tip is free here
    np.testing.assert_allclose(forces[1], [10 * 0.05, 0.0, 0.0])  # a wide tip is pushed out
    # Inside the capillary radius both are pushed, each by its own penetration
    forces = exclusion.calculate_forces_vectorized(tips([0.05, 0.05], [0.002, 0.005]))
    np.testing.assert_allclose(forces[:, 0], [10 * 0.05, 10 * 0.15])


def test_wide_radius_zero_means_one_radius_for_every_tip():
    exclusion, _ = build(wide_radius=0.0)
    forces = exclusion.calculate_forces_vectorized(tips([0.15, 0.05], [0.005, 0.005]))
    np.testing.assert_allclose(forces[0], 0.0)
    np.testing.assert_allclose(forces[1], [10 * 0.05, 0.0, 0.0])


def test_no_demand_site_lies_inside_the_fovea():
    _, demand = build()
    distance = np.hypot(demand.sites[:, 0], demand.sites[:, 1])
    assert distance.min() > 0.1
    # ... and the sites just outside are still there (the lattice was only clipped)
    assert ((distance > 0.1) & (distance < 0.2)).sum() > 0

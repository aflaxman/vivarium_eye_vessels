"""Tests for the ArcadeGuidance astrocyte-template force."""

import copy

import numpy as np
import pandas as pd
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import ArcadeGuidance
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
    "arcade_guidance": {"magnitude": 0.3, "min_radius": 0.006},
}


def build_guidance(**overrides) -> ArcadeGuidance:
    config = copy.deepcopy(CONFIGURATION)
    config["arcade_guidance"].update(overrides)
    guidance = ArcadeGuidance()
    InteractiveContext(
        components=[Particle3D(), PathFreezer(), guidance], configuration=config
    )
    return guidance


def tips(radii) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "x": [1.0, 1.0, 1.5, 1.0],
            "y": [0.0, 0.5, 0.0, -0.3],
            "z": 0.0,
            "frozen": False,
            "path_id": 1,
            "radius": radii,
        }
    )


def test_wide_tips_are_pushed_away_from_the_disc_and_thin_tips_are_not():
    guidance = build_guidance()
    frame = tips([0.02, 0.01, 0.006, 0.003])
    forces = guidance.calculate_forces_vectorized(frame)
    # At the disc itself there is no direction to push
    np.testing.assert_allclose(forces[0], 0.0)
    # Straight up from the disc, straight right of it: unit radial x magnitude
    np.testing.assert_allclose(forces[1], [0.0, 0.3, 0.0], atol=1e-12)
    np.testing.assert_allclose(forces[2], [0.3, 0.0, 0.0], atol=1e-12)
    # A capillary tip is left to the hypoxia signal
    np.testing.assert_allclose(forces[3], 0.0)


def test_zero_magnitude_disables_the_guidance():
    guidance = build_guidance(magnitude=0.0)
    np.testing.assert_allclose(guidance.calculate_forces_vectorized(tips([0.02] * 4)), 0.0)

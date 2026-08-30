"""Tests for Murray's-law vessel calibers (roadmap idea 1)."""

import numpy as np
import pandas as pd
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import EllipsoidContainment
from vivarium_eye_vessels.components.particles import (
    Particle3D,
    PathFreezer,
    PathSplitter,
    murray_bifurcation_angles,
    murray_daughter_radii,
)
from vivarium_eye_vessels.vnv import metrics

CONFIGURATION = {
    "randomness": {"random_seed": 42},
    "time": {
        "start": {"year": 2025, "month": 1, "day": 1},
        "end": {"year": 2025, "month": 3, "day": 1},
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
        "max_depth": 3,
        "murray_exponent": 3.0,
        "flow_asymmetry": 0.15,
        "min_radius": 0.002,
    },
    "ellipsoid_containment": {
        "a": 2,
        "b": 2,
        "c": 0.2,
        "force_type": "hookean",
        "spring_constant": 3,
    },
}


def test_murray_daughter_radii_satisfy_murrays_law():
    major, minor = murray_daughter_radii(1.0, 0.4, 3.0)
    np.testing.assert_allclose(major**3 + minor**3, 1.0, atol=1e-12)
    assert minor < major < 1.0


def test_murray_angles_symmetric_case():
    # Equal daughters under k=3: each deviates ~37.5 degrees, ~75 degrees total
    r_daughter = 2 ** (-1 / 3)
    theta_1, theta_2 = murray_bifurcation_angles(1.0, r_daughter, r_daughter)
    np.testing.assert_allclose(np.degrees(theta_1), 37.47, atol=0.1)
    np.testing.assert_allclose(theta_1, theta_2)


def test_murray_angles_asymmetric_case():
    # The small daughter branches at a wider angle than the large one
    major, minor = murray_daughter_radii(1.0, 0.2, 3.0)
    theta_major, theta_minor = murray_bifurcation_angles(1.0, major, minor)
    assert theta_major < theta_minor


def test_junction_exponents_recover_murray_exponent():
    major, minor = murray_daughter_radii(0.02, 0.4, 3.0)
    pop = pd.DataFrame(
        {
            "x": [0.0, 1.0, 0.0],
            "y": [0.0, 0.0, 1.0],
            "z": [0.0, 0.0, 0.0],
            "frozen": [True, False, False],
            "path_id": [0, 1, 1],
            "parent_id": [-1, 0, 0],
            "depth": [0, 1, 1],
            "radius": [0.02, major, minor],
        },
        index=[0, 1, 2],
    )
    exponents = metrics.junction_exponents(pop)
    np.testing.assert_allclose(exponents, [3.0], atol=1e-6)


def test_simulated_network_has_murray_calibers():
    sim = InteractiveContext(
        components=[Particle3D(), PathFreezer(), PathSplitter(), EllipsoidContainment()],
        configuration=CONFIGURATION,
    )
    for _ in range(150):
        sim.step()
    pop = sim.get_population(["x", "y", "z", "frozen", "path_id", "parent_id", "radius"])

    on_path = pop[pop.path_id >= 0]
    assert (on_path.radius > 0).all(), "vessel particles should carry a caliber"
    assert on_path.radius.max() <= 0.02 + 1e-12, "calibers should not exceed the root radius"

    # Calibers never grow along tree edges (taper <= 1 and Murray daughters shrink),
    # except where the caliber floor binds at min_radius
    children = on_path[on_path.parent_id >= 0]
    children = children[children.parent_id.isin(pop.index)]
    parent_radii = pop.loc[children.parent_id, "radius"].values
    min_radius = CONFIGURATION["path_splitter"]["min_radius"]
    grew = children.radius.values > parent_radii + 1e-12
    at_floor = children.radius.values <= min_radius + 1e-12
    assert not np.any(grew & ~at_floor), "calibers grew along a tree edge"

    # Bifurcations should carry the configured Murray exponent
    exponents = metrics.junction_exponents(pop)
    assert len(exponents) > 0, "no calibered bifurcations found"
    np.testing.assert_allclose(np.median(exponents), 3.0, atol=0.1)

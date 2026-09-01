"""Tests for comb-like side branching off wide trunks and its V&V metric."""

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
from vivarium_eye_vessels.vnv import metrics

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
        "caliber_cadence_exponent": 1.0,
        "side_branch_flow": 0.1,
        "side_branch_radius": 0.008,
        "side_branch_probability": 0.9,
    },
    "ellipsoid_containment": {
        "a": 2,
        "b": 2,
        "c": 0.2,
        "force_type": "hookean",
        "spring_constant": 3,
    },
}


def make_splitter() -> PathSplitter:
    splitter = PathSplitter()
    InteractiveContext(
        components=[Particle3D(), PathFreezer(), splitter, EllipsoidContainment()],
        configuration=copy.deepcopy(CONFIGURATION),
    )
    return splitter


def test_side_branching_trunks_are_exempt_from_cadence_damping():
    splitter = make_splitter()
    active = pd.DataFrame({"radius": [0.02, 0.004, 0.002]}, index=[1, 2, 3])
    probabilities = splitter.split_probabilities(active)
    # The wide trunk side-branches at the full cadence despite exponent 1.0
    np.testing.assert_allclose(probabilities[1], 0.9)
    # Below side_branch_radius the caliber cadence still applies
    np.testing.assert_allclose(probabilities[2], 0.45)
    np.testing.assert_allclose(probabilities[3], 0.9)


def test_side_branch_split_is_asymmetric_and_near_perpendicular():
    splitter = make_splitter()
    pop = pd.DataFrame({"radius": [0.02, 0.003]}, index=[7, 11])
    major, minor, angle_major, angle_minor = splitter.split_radii_and_angles(
        pop, pd.Index([7, 11])
    )

    # Wide trunk: the major daughter continues at nearly the parent caliber,
    # the tooth takes the Murray caliber for a flow fraction of ~0.075-0.125
    assert 0.95 < major[7] / 0.02 < 0.98
    assert 0.40 < minor[7] / 0.02 < 0.52
    # The tooth leaves much more steeply than the trunk deviates
    assert abs(angle_minor[7]) > np.radians(35)
    assert abs(angle_minor[7]) > 2 * abs(angle_major[7])

    # Capillary parent: ordinary near-symmetric dichotomy
    assert 0.7 < minor[11] / major[11] <= 1.0


def test_side_branching_waits_for_its_start_time():
    """Before side_branch_start_time the splitter behaves caliber-blind."""
    config = copy.deepcopy(CONFIGURATION)
    config["path_splitter"]["side_branch_start_time"] = "2025-03-01"  # after sim start
    splitter = PathSplitter()
    InteractiveContext(
        components=[Particle3D(), PathFreezer(), splitter, EllipsoidContainment()],
        configuration=config,
    )
    active = pd.DataFrame({"radius": [0.02]}, index=[1])
    # The cadence damping applies to the wide tip: no comb exemption yet
    np.testing.assert_allclose(splitter.split_probabilities(active)[1], 0.09)
    # And splits are ordinary near-symmetric dichotomies
    pop = pd.DataFrame({"radius": [0.02]}, index=[7])
    major, minor, _, _ = splitter.split_radii_and_angles(pop, pd.Index([7]))
    assert minor[7] / major[7] > 0.7


def test_anastomosis_targets_respect_min_layer():
    tips = pd.DataFrame({"x": [0.0], "y": [0.0], "z": [0.0], "vessel_type": [1]}, index=[0])
    frozen = pd.DataFrame(
        {
            "x": [0.01, 0.02],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0],
            "vessel_type": [2, 2],
            "radius": [0.003, 0.003],
            "layer_id": [0, 1],
        },
        index=[10, 11],
    )
    from vivarium_eye_vessels.components.particles import anastomosis_targets

    neighbor_lists = [[0, 1]]
    # Without the layer floor the nearer (superficial) target wins
    all_layers = anastomosis_targets(tips, frozen, neighbor_lists, 0.004)
    assert all_layers[0] == 10
    # With min_layer 1 only the deep target qualifies
    deep_only = anastomosis_targets(tips, frozen, neighbor_lists, 0.004, min_layer=1)
    assert deep_only[0] == 11


def test_wide_junction_spacing_counts_comb_teeth():
    """A synthetic comb: thick trunk with teeth every 20 px reads ~20 px spacing."""
    image = np.zeros((60, 220), dtype=bool)
    image[26:35, 10:210] = True  # 9 px wide trunk -> diameter > 4 px
    for x in range(20, 200, 20):
        image[35:55, x : x + 2] = True  # thin teeth, diameter ~2 px

    from skimage.morphology import skeletonize

    spacing = metrics.wide_junction_spacing(skeletonize(image), image)
    assert 12 < spacing < 28

    # Without teeth there are no junctions: the whole trunk is one span
    bare = np.zeros((60, 220), dtype=bool)
    bare[26:35, 10:210] = True
    bare_spacing = metrics.wide_junction_spacing(skeletonize(bare), bare)
    assert bare_spacing > 150

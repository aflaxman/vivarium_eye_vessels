"""Unit tests for the V&V network metrics (no downloads, no display)."""

import numpy as np
import pandas as pd

from vivarium_eye_vessels.vnv import metrics


def test_box_counting_dimension_of_line_is_near_one():
    image = np.zeros((256, 256), dtype=bool)
    image[128, :] = True
    dimension = metrics.box_counting_dimension(image)
    assert 0.9 < dimension < 1.1


def test_box_counting_dimension_of_filled_square_is_near_two():
    image = np.ones((256, 256), dtype=bool)
    dimension = metrics.box_counting_dimension(image)
    assert 1.9 < dimension < 2.1


def test_skeleton_branches_straight_line():
    image = np.zeros((64, 64), dtype=bool)
    image[32, 8:56] = True
    branches = metrics.skeleton_branches(image)
    assert len(branches) == 1
    np.testing.assert_allclose(branches.tortuosity.iloc[0], 1.0, atol=0.05)


def test_rasterize_network_draws_segments():
    edges = pd.DataFrame(
        {
            "x0": [-1.0],
            "y0": [0.0],
            "z0": [0.0],
            "x1": [1.0],
            "y1": [0.0],
            "z1": [0.0],
            "child": [1],
            "parent": [0],
        }
    )
    image = metrics.rasterize_network(edges, bounds=(2.0, 2.0), size=128)
    assert image.sum() > 50  # a horizontal line across half the image


def test_bifurcation_angles_right_angle():
    pop = pd.DataFrame(
        {
            "x": [0.0, 1.0, 0.0],
            "y": [0.0, 0.0, 1.0],
            "z": [0.0, 0.0, 0.0],
            "frozen": [True, False, False],
            "path_id": [0, 1, 1],
            "parent_id": [-1, 0, 0],
            "depth": [0, 1, 1],
        },
        index=[0, 1, 2],
    )
    angles = metrics.bifurcation_angles(pop)
    np.testing.assert_allclose(angles, [90.0], atol=1e-6)


def test_path_tortuosity_straight_chain_is_one():
    pop = pd.DataFrame(
        {
            "x": [0.0, 1.0, 2.0, 3.0],
            "y": [0.0, 0.0, 0.0, 0.0],
            "z": [0.0, 0.0, 0.0, 0.0],
            "frozen": [True, True, True, False],
            "path_id": [7, 7, 7, 7],
            "parent_id": [-1, 0, 1, 2],
            "depth": [0, 0, 0, 0],
        },
        index=[0, 1, 2, 3],
    )
    ratios = metrics.path_tortuosity(pop)
    np.testing.assert_allclose(ratios, [1.0], atol=1e-9)

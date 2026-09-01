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


def test_skeleton_branch_diameter_recovers_bar_width():
    from skimage.morphology import skeletonize

    image = np.zeros((64, 128), dtype=bool)
    image[28:37, 8:120] = True  # a 9-px-thick horizontal bar
    branches = metrics.skeleton_branches(skeletonize(image), image)
    assert len(branches) == 1
    assert 7.0 < branches.diameter_px.iloc[0] < 11.0


def test_stratify_by_diameter_bins():
    from vivarium_eye_vessels.vnv.compare import stratify_by_diameter

    # A 1-px skeleton line measures diameter exactly 2.0, so the capillary
    # bin is closed at 2 px
    strata = stratify_by_diameter([10.0, 20.0, 30.0, 40.0], [1.0, 2.0, 3.9, 5.0])
    np.testing.assert_allclose(strata["diameter_le_2px"], [10.0, 20.0])
    np.testing.assert_allclose(strata["diameter_2_4px"], [30.0])
    np.testing.assert_allclose(strata["diameter_gt_4px"], [40.0])


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


def test_skeleton_pixel_diameters_recover_bar_widths():
    """A wide bar and a thin bar contribute pixels at their own diameters."""
    image = np.zeros((40, 220), dtype=bool)
    image[15:24, 10:110] = True  # 9 px wide bar
    image[30:32, 10:210] = True  # 2 px thin bar
    diameters = metrics.skeleton_pixel_diameters(image)
    assert len(diameters) > 100
    wide = diameters[diameters > 6]
    thin = diameters[diameters <= 4]
    # The thin bar is twice as long, so it contributes more skeleton pixels
    assert len(thin) > len(wide) > 30
    assert 7.0 < np.median(wide) < 11.0
    assert 1.5 < np.median(thin) < 3.5

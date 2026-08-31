"""Tests for the calibration objective and search (roadmap idea 8)."""

import numpy as np

from vivarium_eye_vessels.vnv.calibrate import (
    SEARCH_SPACE,
    TARGETS,
    calibration_score,
    combine_seed_scores,
    coordinate_descent,
)


def on_target_stats() -> dict:
    return {name: spec["target"] for name, spec in TARGETS.items()}


def test_on_target_stats_score_zero():
    scores = calibration_score(on_target_stats())
    assert scores["total"] == 0.0


def test_deviation_scores_as_squared_z():
    stats = on_target_stats()
    stats["skeleton_density"] = (
        TARGETS["skeleton_density"]["target"] + 2 * TARGETS["skeleton_density"]["scale"]
    )
    scores = calibration_score(stats)
    np.testing.assert_allclose(scores["skeleton_density"], 4.0)
    np.testing.assert_allclose(scores["total"], 4.0)


def test_one_sided_targets_only_penalize_the_bad_direction():
    stats = on_target_stats()
    # Better-than-target perfusion and KS cost nothing
    stats["perfused_fraction"] = 1.0
    stats["ks_log_length"] = 0.0
    scores = calibration_score(stats)
    assert scores["perfused_fraction"] == 0.0
    assert scores["ks_log_length"] == 0.0
    # Worse-than-target perfusion is penalized
    stats["perfused_fraction"] = 0.94
    np.testing.assert_allclose(calibration_score(stats)["perfused_fraction"], 4.0)


def test_missing_or_nan_stats_never_win():
    stats = on_target_stats()
    stats["artery_vein_caliber_ratio"] = float("nan")
    del stats["fractal_dimension"]
    scores = calibration_score(stats)
    assert scores["artery_vein_caliber_ratio"] == 25.0
    assert scores["fractal_dimension"] == 25.0


def test_multiseed_objective_is_the_mean_across_seeds():
    per_seed = {
        7: {"skeleton_density": 1.0, "wide_share": 3.0, "total": 4.0},
        42: {"skeleton_density": 5.0, "wide_share": 1.0, "total": 6.0},
    }
    combined = combine_seed_scores(per_seed)
    np.testing.assert_allclose(combined["skeleton_density"], 3.0)
    np.testing.assert_allclose(combined["wide_share"], 2.0)
    np.testing.assert_allclose(combined["total"], 5.0)


def test_multiseed_objective_punishes_a_single_collapsed_seed():
    steady = combine_seed_scores({1: {"total": 60.0}, 2: {"total": 60.0}})
    lucky_but_fragile = combine_seed_scores({1: {"total": 40.0}, 2: {"total": 1800.0}})
    assert steady["total"] < lucky_but_fragile["total"]


def test_coordinate_descent_finds_the_known_optimum():
    space = {("a", "x"): [1, 2, 3], ("b", "y"): [10, 20, 30]}
    optimum = {("a", "x"): 2, ("b", "y"): 30}
    calls = []

    def evaluate(overrides):
        calls.append(dict(overrides))
        return sum(
            (overrides.get(knob, values[0]) - optimum[knob]) ** 2
            for knob, values in space.items()
        )

    outcome = coordinate_descent(evaluate, space, {}, budget=30)
    assert outcome["best"][("a", "x")] == 2
    assert outcome["best"][("b", "y")] == 30
    assert outcome["best_score"] == 0
    assert outcome["evaluations"] <= 30


def test_coordinate_descent_respects_budget():
    def evaluate(overrides):
        return 1.0  # nothing ever improves

    outcome = coordinate_descent(evaluate, SEARCH_SPACE, {}, budget=5)
    assert outcome["evaluations"] <= 5

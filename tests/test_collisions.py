"""Tips do not pass through vessels of their own plane (CollisionGuard)."""

from types import SimpleNamespace

import numpy as np
import pandas as pd

from vivarium_eye_vessels.components.collisions import CollisionGuard, segments_cross


def make_guard(**overrides) -> CollisionGuard:
    guard = CollisionGuard()
    settings = dict(
        enabled=True,
        layers=[0],
        capillary_radius=0.00095,
        crossing_min_radius=0.005,
        search_radius=0.05,
        mode="slide",
        crossing_wide_ratio=0.0,
    )
    settings.update(overrides)
    guard.config = SimpleNamespace(**settings)
    guard.enabled = True
    guard.layers = settings["layers"]
    guard.capillary_radius = settings["capillary_radius"]
    guard.crossing_min_radius = settings["crossing_min_radius"]
    guard.search_radius = settings["search_radius"]
    guard.mode = settings["mode"]
    guard.crossing_wide_ratio = settings["crossing_wide_ratio"]
    guard.total_blocked = 0
    return guard


def test_segments_cross_is_a_proper_intersection():
    p0, p1 = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]), np.array(
        [[1.0, 1.0], [1.0, 1.0], [1.0, 0.0]]
    )
    q0, q1 = np.array([[0.0, 1.0], [2.0, 0.0], [0.5, 0.0]]), np.array(
        [[1.0, 0.0], [3.0, 1.0], [0.5, 1.0]]
    )
    # crossing X; disjoint; T touching at an endpoint of q (u = 0 -> not proper)
    assert list(segments_cross(p0, p1, q0, q1)) == [True, False, False]


def population():
    # a frozen horizontal vessel (0 -> 1 -> 2) of tree 1 along y = 0, a frozen vessel of tree 2
    # (3 -> 4) at y = 0.1, two frozen stubs (5, 6) the tips descend from, and tips 10-14
    return pd.DataFrame(
        {
            "x": [0.0, 0.02, 0.04, 0.0, 0.04, 0.01, 0.02, 0.01, 0.01, 0.02, 0.02, 0.035],
            "y": [0.0, 0.0, 0.0, 0.1, 0.1, -0.05, 0.05, -0.005, -0.005, 0.095, 0.095, -0.005],
            "frozen": [True] * 7 + [False] * 5,
            "parent_id": [-1, 0, 1, -1, 3, -1, -1, 5, 5, 6, 6, 2],
            "path_id": [0, 0, 0, 1, 1, 5, 6, 7, 8, 9, 9, 11],
            "radius": [
                0.006,
                0.006,
                0.006,
                0.006,
                0.006,
                0.006,
                0.006,
                0.003,
                0.003,
                0.006,
                0.003,
                0.003,
            ],
            "vessel_type": [1, 1, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1],
            "layer_id": 0,
        },
        index=[0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14],
    )


def test_same_tree_crossings_are_blocked_and_wide_opposite_trees_pass():
    pop = population()
    guard = make_guard()
    tips = pop.loc[[10, 11, 12, 13, 14]]
    after = tips[["x", "y"]].copy()
    after["y"] = [
        0.005,
        0.005,
        0.105,
        0.105,
        0.005,
    ]  # every tip steps across the vessel above it
    blocked = set(guard.blocked_tips(tips, after, pop))
    # 10: tree 1 crossing tree 1 -> blocked; 11: tree 2 (3 px) crossing tree 1 -> too thin -> blocked
    # 12: tree 1 (6 px) crossing tree 2 (6 px) -> both wide -> allowed
    # 13: tree 2 crossing tree 2 -> blocked; 14: its parent is particle 2 and the segment it crosses
    # is 1 -> 2, its own trail -> not an obstacle
    assert blocked == {10, 11, 13}


def test_guard_ignores_other_layers_and_capillaries():
    pop = population()
    pop.loc[[10], "layer_id"] = 2
    pop.loc[[11], "radius"] = 0.0009
    guard = make_guard()
    tips = pop.loc[[10, 11]]
    after = tips[["x", "y"]].copy()
    after["y"] = 0.005
    assert len(guard.blocked_tips(tips, after, pop)) == 0


def test_a_blocked_tip_slides_along_the_vessel_or_stops():
    pop = population()
    # tip 10 heads up and to the right into the horizontal vessel at y = 0
    updates = pd.DataFrame(
        {"x": [0.02], "y": [0.005], "z": [0.0], "vx": [0.2], "vy": [0.2], "vz": [0.0]},
        index=[10],
    )
    before = pop.loc[[10], ["x", "y"]].assign(z=0.0)
    guard = make_guard(mode="slide")
    blocked = guard.deflect(before, updates, pop)
    assert list(blocked) == [10]
    # only the component along the vessel (x) survives; it does not cross y = 0
    assert np.isclose(updates.loc[10, "y"], -0.005) and updates.loc[10, "x"] > 0.01
    assert np.isclose(updates.loc[10, "vy"], 0.0) and updates.loc[10, "vx"] > 0
    updates = pd.DataFrame(
        {"x": [0.02], "y": [0.005], "z": [0.0], "vx": [0.2], "vy": [0.2], "vz": [0.0]},
        index=[10],
    )
    guard = make_guard(mode="stop")
    guard.deflect(before, updates, pop)
    assert np.isclose(updates.loc[10, "x"], 0.01) and np.isclose(updates.loc[10, "y"], -0.005)


def test_a_blocked_tip_can_reflect_off_the_vessel():
    pop = population()
    updates = pd.DataFrame(
        {"x": [0.02], "y": [0.005], "z": [0.0], "vx": [0.2], "vy": [0.2], "vz": [0.0]},
        index=[10],
    )
    before = pop.loc[[10], ["x", "y"]].assign(z=0.0)
    guard = make_guard(mode="reflect")
    assert list(guard.deflect(before, updates, pop)) == [10]
    # the move along the vessel is kept, the move toward it reversed: (0.01, -0.005) + (0.01, -0.01)
    assert np.isclose(updates.loc[10, "x"], 0.02) and np.isclose(updates.loc[10, "y"], -0.015)
    assert np.isclose(updates.loc[10, "vx"], 0.2) and np.isclose(updates.loc[10, "vy"], -0.2)


def test_a_much_wider_tip_passes_over_a_twig():
    pop = population()
    pop.loc[[0, 1, 2], "radius"] = 0.002  # the tree-1 vessel is a twig
    pop.loc[10, "radius"] = 0.006  # a trunk of the same tree heads across it
    tips = pop.loc[[10]]
    after = tips[["x", "y"]].copy()
    after["y"] = 0.005
    assert list(make_guard().blocked_tips(tips, after, pop)) == [10]
    assert len(make_guard(crossing_wide_ratio=2.0).blocked_tips(tips, after, pop)) == 0

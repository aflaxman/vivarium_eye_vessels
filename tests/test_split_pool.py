"""Tests for the PathSplitter free-pool carving, sprout path ids, in-plane
split geometry, and OU-state hand-off — the defects the bug hunt confirmed."""

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
from vivarium_eye_vessels.vnv import simulation

CONFIGURATION = {
    "randomness": {"random_seed": 7},
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
        "split_interval": 10,
        "split_angle": 60,
        "split_probability": 0.9,
        "max_depth": 6,
        "murray_exponent": 3.0,
        "flow_asymmetry": 0.15,
        "min_radius": 0.002,
        # Force the re-sprout phase to fire alongside the active-split phase,
        # so both draw from the free pool in the same round
        "min_active_tips": 5,
        "resprout_established_size": 1,
    },
    "ellipsoid_containment": {
        "a": 2,
        "b": 2,
        "c": 0.2,
        "force_type": "hookean",
        "spring_constant": 3,
    },
}


def build(**splitter_overrides):
    config = copy.deepcopy(CONFIGURATION)
    config["path_splitter"].update(splitter_overrides)
    splitter = PathSplitter()
    sim = InteractiveContext(
        components=[Particle3D(), PathFreezer(), splitter, EllipsoidContainment()],
        configuration=config,
    )
    return sim, splitter


def test_split_phases_never_claim_a_free_slot_twice():
    sim, splitter = build()
    seen = []
    original_commit = splitter.commit

    def spy(updates):
        if updates:
            combined = pd.concat(updates, axis=0)
            seen.append(bool(combined.index.duplicated().any()))
        return original_commit(updates)

    splitter.commit = spy
    simulation.run_steps(sim, 120)
    assert seen, "no split rounds ran"
    assert not any(seen), "a split round wrote two branches to the same particle"


def test_sprouts_get_distinct_path_ids():
    sim, splitter = build()
    pop = sim.get_population(splitter.required_attributes)
    frozen = pop[pop.frozen & (pop.path_id >= 0)]
    to_split = frozen.index[:4]
    available = pop[~pop.frozen & (pop.path_id < 0)]
    before = splitter.next_path_id
    updates = splitter.split_frozen(pop, to_split, available) or []
    sprout_paths = [f.path_id.iloc[0] for f in updates]
    assert len(sprout_paths) == len(set(sprout_paths)), "sprouts shared a path id"
    assert splitter.next_path_id == before + len(updates), "next_path_id not advanced"


def test_split_stays_in_the_plexus_plane():
    _, splitter = build()
    # A tip heading almost along +x with a small out-of-plane component: the
    # old axis threw the daughters into z; the fix keeps the split horizontal
    for vel in ([0.3, 0.02, 0.02], [0.3, 0.0, 0.0], [0.0, 0.3, 0.01]):
        vel = np.array(vel, dtype=float)
        axis = splitter.split_axis(vel / np.linalg.norm(vel))
        rot = splitter._rotation_matrix(axis, np.radians(30))
        daughter = rot @ vel
        # the daughter's z stays within the parent's z scale, not amplified by
        # rotating a large in-plane speed out of plane
        assert abs(daughter[2]) <= abs(vel[2]) + 1e-9


def test_freezer_carries_the_ou_steering_state():
    sim, _ = build()
    simulation.run_steps(sim, 6)
    components = sim.list_components()
    freezer = next(c for c in components.values() if isinstance(c, PathFreezer))
    pop = sim.get_population(freezer.required_attributes)
    tip = pop[~pop.frozen & (pop.path_id >= 0)].index[0]
    pop.loc[tip, ["wx", "wy", "wz"]] = [0.11, -0.22, 0.33]
    freezer.freeze_particles(pop)  # hands the tip's identity to a free particle
    after = sim.get_population(["wx", "wy", "wz", "frozen", "path_id", "parent_id"])
    continuation = after[(after.parent_id == tip) & ~after.frozen & (after.path_id >= 0)]
    assert not continuation.empty, "the frozen tip produced no continuation"
    row = continuation.iloc[0]
    np.testing.assert_allclose([row.wx, row.wy, row.wz], [0.11, -0.22, 0.33])

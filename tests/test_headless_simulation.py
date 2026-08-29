"""Headless integration tests for the vessel growth components.

These tests exercise the particle system, path freezing, splitting, and
boundary forces without the pygame visualizer, so they run on CI machines
with no display.
"""

import numpy as np
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import (
    EllipsoidContainment,
    FrozenRepulsion,
)
from vivarium_eye_vessels.components.particles import (
    Particle3D,
    PathExtinction,
    PathFreezer,
    PathSplitter,
)

CONFIGURATION = {
    "randomness": {"random_seed": 12345},
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
    },
    "path_freezer": {"freeze_interval": 3},
    "path_splitter": {
        "split_interval": 20,
        "split_angle": 60,
        "split_probability": 0.9,
        "max_depth": 2,
    },
    "path_extinction": {"force_threshold": 0.85},
    "frozen_repulsion": {
        "force_type": "hookean",
        "spring_constant": 1.5,
        "interaction_radius": 0.15,
        "delay": 1,
    },
    "ellipsoid_containment": {
        "a": 2,
        "b": 2,
        "c": 0.2,
        "force_type": "hookean",
        "spring_constant": 3,
    },
}

PARTICLE_ATTRIBUTES = [
    "x",
    "y",
    "z",
    "vx",
    "vy",
    "vz",
    "frozen",
    "path_id",
    "parent_id",
    "depth",
]


def make_simulation() -> InteractiveContext:
    return InteractiveContext(
        components=[
            Particle3D(),
            PathFreezer(),
            PathExtinction(),
            PathSplitter(),
            EllipsoidContainment(),
            FrozenRepulsion(),
        ],
        configuration=CONFIGURATION,
    )


def test_simulation_setup_and_initialization():
    sim = make_simulation()
    pop = sim.get_population(PARTICLE_ATTRIBUTES)

    assert len(pop) == 300
    assert (pop.path_id >= 0).sum() == 4, "expected 4 initial vessel tips"
    assert not pop.frozen.any(), "no particles should start frozen"


def test_vessels_grow_and_branch():
    sim = make_simulation()
    for _ in range(100):
        sim.step()
    pop = sim.get_population(PARTICLE_ATTRIBUTES)

    frozen = pop[pop.frozen]
    assert len(frozen) > 0, "no vessel segments were created"

    on_path = pop[pop.path_id >= 0]
    assert len(on_path) > 4, "vessel paths did not grow"

    has_parents = pop[pop.parent_id >= 0]
    assert len(has_parents) > 0, "no parent-child relationships formed"

    # PathSplitter should have created deeper branches by now
    assert pop.depth.max() >= 1, "no branching occurred"
    assert pop.depth.max() <= 2, "max_depth exceeded"


def test_ellipsoid_containment_bounds_particles():
    sim = make_simulation()
    for _ in range(100):
        sim.step()
    pop = sim.get_population(["x", "y", "z"])

    ellipsoid_values = (
        (pop.x / 2.0) ** 2 + (pop.y / 2.0) ** 2 + (pop.z / 0.2) ** 2
    ).to_numpy()
    # Hookean containment is soft, so allow modest excursions past the surface
    assert np.quantile(ellipsoid_values, 0.99) < 2.0, "particles escaped the ellipsoid"


def test_more_particles_are_added_when_needed():
    sim = make_simulation()
    for _ in range(250):
        sim.step()
    pop = sim.get_population(["frozen"])

    assert len(pop) > 300, "simulant creator never added new particles"

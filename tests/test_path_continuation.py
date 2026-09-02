"""Tests for continuation_keeps_path: the continuing daughter of a split keeps
its parent's path_id, so its own fresh trail stays under FrozenRepulsion's
same-path delay exemption."""

import copy

from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import EllipsoidContainment
from vivarium_eye_vessels.components.particles import (
    Particle3D,
    PathFreezer,
    PathSplitter,
)
from vivarium_eye_vessels.vnv import simulation

N_VESSELS = 4

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
        "initial_circle": {"center": [1.0, 0.0, 0.0], "radius": 0.05, "n_vessels": N_VESSELS},
        "root_radius": 0.02,
    },
    "path_freezer": {"freeze_interval": 3, "radius_taper": 0.999},
    "path_splitter": {
        "split_interval": 20,
        "split_angle": 60,
        "split_probability": 1.0,
        "max_depth": 6,
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


def depth_zero_paths(keeps_path: bool) -> set[int]:
    config = copy.deepcopy(CONFIGURATION)
    config["path_splitter"]["continuation_keeps_path"] = keeps_path
    sim = InteractiveContext(
        components=[Particle3D(), PathFreezer(), PathSplitter(), EllipsoidContainment()],
        configuration=config,
    )
    simulation.run_steps(sim, 45)  # past the split rounds at steps 20 and 40
    pop = sim.get_population(["frozen", "path_id", "depth"])
    trunks = pop[~pop.frozen & (pop.path_id >= 0) & (pop.depth == 0)]
    assert not trunks.empty
    return set(trunks.path_id)


def test_continuing_daughter_keeps_the_root_path():
    # Every surviving trunk tip still carries one of the root path ids
    assert depth_zero_paths(keeps_path=True) <= set(range(N_VESSELS))


def test_legacy_relabels_the_continuing_daughter():
    # Legacy: after its first split a trunk's continuation is a new path
    assert depth_zero_paths(keeps_path=False) - set(range(N_VESSELS))


def test_freezer_never_ends_the_front_when_the_pool_runs_dry():
    # A pool of 20 with 4 roots runs dry after a few freeze rounds; the
    # freezer must top it up and skip the round rather than freeze every tip
    # without a continuation (which silently ended the whole front)
    config = copy.deepcopy(CONFIGURATION)
    config["population"]["population_size"] = 20
    sim = InteractiveContext(
        components=[Particle3D(), PathFreezer(), PathSplitter(), EllipsoidContainment()],
        configuration=config,
    )
    for _ in range(30):
        sim.step()
        pop = sim.get_population(["frozen", "path_id"])
        assert (~pop.frozen & (pop.path_id >= 0)).sum() >= N_VESSELS

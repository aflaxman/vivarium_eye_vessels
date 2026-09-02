"""Tests for the DevelopmentalWave closed-loop growth front."""

import copy
from types import SimpleNamespace

import numpy as np
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import (
    DevelopmentalWave,
    EllipsoidContainment,
    PerfusionDemand,
)
from vivarium_eye_vessels.components.particles import (
    VESSEL_TYPE_ARTERY,
    Particle3D,
    PathFreezer,
    PathSplitter,
)
from vivarium_eye_vessels.vnv import simulation

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
        "max_depth": 6,
        "murray_exponent": 3.0,
        "flow_asymmetry": 0.15,
        "min_radius": 0.002,
    },
    "perfusion_demand": {
        "site_spacing": 0.25,
        "perfusion_radius": 0.15,
        "influence_radius": 2.0,
        "magnitude": 0.3,
    },
    "developmental_wave": {
        "enabled": True,
        "start_radius": 0.1,
        "wave_speed": 0.01,
        "lookahead": 0.3,
        "advance_threshold": 0.5,
        "hold_resprout_steps": 2,
        "resprout_count": 2,
    },
    "ellipsoid_containment": {
        "a": 2,
        "b": 2,
        "c": 0.2,
        "force_type": "hookean",
        "spring_constant": 3,
    },
}


def build_sim(**overrides) -> tuple[InteractiveContext, DevelopmentalWave]:
    config = copy.deepcopy(CONFIGURATION)
    for key, value in overrides.items():
        config["developmental_wave"][key] = value
    wave = DevelopmentalWave()
    sim = InteractiveContext(
        components=[
            Particle3D(),
            PathFreezer(),
            PathSplitter(),
            EllipsoidContainment(),
            PerfusionDemand(),
            wave,
        ],
        configuration=config,
    )
    return sim, wave


def test_visible_exposes_only_sites_near_the_front():
    _, wave = build_sim()
    sites = np.array([[1.0, 0.0, 0.0], [1.35, 0.0, 0.0], [1.45, 0.0, 0.0]])
    # radius 0.1 + lookahead 0.3 admits the first two, not the third
    np.testing.assert_array_equal(wave.visible(sites), sites[:2])

    _, disabled = build_sim(enabled=False)
    np.testing.assert_array_equal(disabled.visible(sites), sites)


def test_wave_advances_when_served_and_holds_when_not():
    sim, wave = build_sim()
    start = wave.radius
    simulation.run_steps(sim, 30)
    assert wave.radius > start, "a served front should advance"

    # Push the front far beyond the network: the tissue inside is unserved,
    # so the front must hold
    wave.radius = 3.0
    simulation.run_steps(sim, 5)
    assert wave.radius == 3.0, "an unserved front must hold"


def test_held_front_resprouts_beside_the_stall():
    sim, wave = build_sim()
    simulation.run_steps(sim, 30)
    pop = sim.get_population(["frozen", "path_id", "vessel_type", "x", "y", "z"])
    frozen = pop[pop.frozen & (pop.path_id >= 0) & (pop.vessel_type == VESSEL_TYPE_ARTERY)]
    assert not frozen.empty

    active_before = (
        ~pop.frozen & (pop.path_id >= 0) & (pop.vessel_type == VESSEL_TYPE_ARTERY)
    ).sum()
    wave.radius = 3.0  # everything inside counts, so the artery tree is stalled
    wave.resprout_toward_stall(VESSEL_TYPE_ARTERY, SimpleNamespace(index=pop.index))

    after = sim.get_population(["frozen", "path_id", "vessel_type", "parent_id"])
    sprouts = after[
        ~after.frozen
        & (after.path_id >= 0)
        & (after.vessel_type == VESSEL_TYPE_ARTERY)
        & after.parent_id.isin(frozen.index)
    ]
    assert len(sprouts) > 0, "a stalled tree should sprout off its frozen frontier"
    total_active = (
        ~after.frozen & (after.path_id >= 0) & (after.vessel_type == VESSEL_TYPE_ARTERY)
    ).sum()
    assert total_active > active_before


def test_disabled_wave_leaves_dynamics_unchanged():
    sim_off, wave_off = build_sim(enabled=False)
    start = wave_off.radius
    simulation.run_steps(sim_off, 45)
    assert wave_off.radius == start

    config = copy.deepcopy(CONFIGURATION)
    del config["developmental_wave"]
    sim_without = InteractiveContext(
        components=[
            Particle3D(),
            PathFreezer(),
            PathSplitter(),
            EllipsoidContainment(),
            PerfusionDemand(),
        ],
        configuration=config,
    )
    simulation.run_steps(sim_without, 45)
    attributes = ["x", "y", "z", "frozen", "path_id", "radius"]
    assert sim_off.get_population(attributes).equals(sim_without.get_population(attributes))

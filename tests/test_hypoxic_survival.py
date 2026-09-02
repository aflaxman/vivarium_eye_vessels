"""Tests for hypoxic survival relief: the extinction-threshold pipeline that
PerfusionDemand raises for tips sitting in unserved tissue."""

import copy

import numpy as np
from scipy.spatial import cKDTree
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import (
    EllipsoidContainment,
    PerfusionDemand,
)
from vivarium_eye_vessels.components.particles import (
    Particle3D,
    PathExtinction,
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
    "path_extinction": {"force_threshold": 2.0},
    "perfusion_demand": {
        "site_spacing": 0.25,
        "perfusion_radius": 0.15,
        "influence_radius": 2.0,
        "magnitude": 0.3,
        "survival_factor": 2.0,
    },
    "ellipsoid_containment": {
        "a": 2,
        "b": 2,
        "c": 0.2,
        "force_type": "hookean",
        "spring_constant": 3,
    },
}


def build_sim(survival_factor: float) -> tuple[InteractiveContext, PerfusionDemand]:
    config = copy.deepcopy(CONFIGURATION)
    config["perfusion_demand"]["survival_factor"] = survival_factor
    demand = PerfusionDemand()
    sim = InteractiveContext(
        components=[
            Particle3D(),
            PathFreezer(),
            PathSplitter(),
            PathExtinction(),
            EllipsoidContainment(),
            demand,
        ],
        configuration=config,
    )
    return sim, demand


def active_tips(sim: InteractiveContext):
    pop = sim.get_population(["x", "y", "z", "frozen", "path_id", "vessel_type"])
    return pop[~pop.frozen & (pop.path_id >= 0)]


def test_relief_raises_the_threshold_exactly_for_tips_in_hypoxic_tissue():
    sim, demand = build_sim(survival_factor=2.0)
    simulation.run_steps(sim, 30)
    tips = active_tips(sim)
    assert not tips.empty
    thresholds = sim.get_value("particle.extinction_threshold")(tips.index)

    # Recompute the predicate independently of the component
    expected = np.zeros(len(tips), dtype=bool)
    for vessel_type in np.unique(tips.vessel_type):
        sites = demand.hypoxic_sites(int(vessel_type))
        selected = (tips.vessel_type == vessel_type).to_numpy()
        distances, _ = cKDTree(sites).query(tips[["x", "y", "z"]].to_numpy()[selected], k=1)
        expected[selected] = distances <= demand.perfusion_radius
    assert expected.any(), "frontier tips should sit in hypoxic tissue"
    np.testing.assert_allclose(thresholds.to_numpy(), np.where(expected, 4.0, 2.0))


def test_factor_one_is_the_legacy_threshold_everywhere():
    sim, _ = build_sim(survival_factor=1.0)
    simulation.run_steps(sim, 30)
    tips = active_tips(sim)
    thresholds = sim.get_value("particle.extinction_threshold")(tips.index)
    assert (thresholds == 2.0).all()

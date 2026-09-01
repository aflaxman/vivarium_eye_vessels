"""Tests for hypoxia-driven growth via space colonization (roadmap idea 2)."""

import copy

import numpy as np
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import (
    EllipsoidContainment,
    FrozenRepulsion,
    PerfusionDemand,
    colonization_forces,
    generate_demand_sites,
)
from vivarium_eye_vessels.components.particles import (
    Particle3D,
    PathExtinction,
    PathFreezer,
    PathSplitter,
)
from vivarium_eye_vessels.vnv import metrics, simulation

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
    },
    "path_extinction": {"force_threshold": 1.2},
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
    "perfusion_demand": {
        "site_spacing": 0.1,
        "perfusion_radius": 0.15,
        "influence_radius": 2.0,
        "magnitude": 0.3,
    },
}

SEMI_AXES = (2.0, 2.0, 0.2)


def test_demand_sites_lie_inside_ellipsoid():
    sites = generate_demand_sites(np.array([2.0, 2.0, 0.2]), 0.1)
    assert len(sites) > 500
    values = (sites[:, 0] / 2) ** 2 + (sites[:, 1] / 2) ** 2 + (sites[:, 2] / 0.2) ** 2
    assert values.max() <= 1.0 + 1e-12


def test_colonization_force_points_toward_hypoxic_tissue():
    tips = np.array([[0.0, 0.0, 0.0]])
    hypoxic = np.array([[1.0, 0.0, 0.0], [1.0, 0.5, 0.0]])
    forces = colonization_forces(tips, hypoxic, influence_radius=2.0, magnitude=0.3)
    np.testing.assert_allclose(np.linalg.norm(forces[0]), 0.3, atol=1e-12)
    assert forces[0, 0] > 0, "force should point toward the hypoxic sites"


def test_each_site_recruits_only_its_nearest_tip():
    tips = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    hypoxic = np.array([[1.9, 0.0, 0.0]])
    forces = colonization_forces(tips, hypoxic, influence_radius=3.0, magnitude=0.3)
    np.testing.assert_allclose(forces[0], 0.0)
    assert forces[1, 0] < 0, "the nearer tip is pulled toward the site"


def test_colonization_force_empty_cases():
    empty = np.zeros((0, 3))
    tips = np.array([[0.0, 0.0, 0.0]])
    assert colonization_forces(empty, tips, 1.0, 0.3).shape == (0, 3)
    np.testing.assert_allclose(colonization_forces(tips, empty, 1.0, 0.3), 0.0)
    # Sites out of influence range exert no pull
    far = np.array([[10.0, 0.0, 0.0]])
    np.testing.assert_allclose(colonization_forces(tips, far, 1.0, 0.3), 0.0)


def run_simulation(with_perfusion: bool, steps: int = 250) -> float:
    """Run the test model and return the perfused tissue fraction."""
    components = [
        Particle3D(),
        PathFreezer(),
        PathExtinction(),
        PathSplitter(),
        EllipsoidContainment(),
        FrozenRepulsion(),
    ]
    if with_perfusion:
        components.append(PerfusionDemand())
    sim = InteractiveContext(
        components=components, configuration=copy.deepcopy(CONFIGURATION)
    )
    simulation.run_steps(sim, steps)
    pop = simulation.get_network(sim)
    return metrics.perfused_fraction(pop, SEMI_AXES, site_spacing=0.1, perfusion_radius=0.15)


def test_caliber_attenuation_scales_the_pull():
    """A wide tip feels (reference/radius)**exponent of the capillary pull."""
    import pandas as pd

    config = copy.deepcopy(CONFIGURATION)
    config["perfusion_demand"]["caliber_exponent"] = 1.0
    demand = PerfusionDemand()
    sim = InteractiveContext(
        components=[
            Particle3D(),
            PathFreezer(),
            PathExtinction(),
            PathSplitter(),
            EllipsoidContainment(),
            FrozenRepulsion(),
            demand,
        ],
        configuration=config,
    )
    simulation.run_steps(sim, 10)

    def pull(radius: float) -> float:
        particles = pd.DataFrame(
            {
                "x": [0.5],
                "y": [0.3],
                "z": [0.0],
                "frozen": [False],
                "path_id": [1],
                "vessel_type": [1],
                "radius": [radius],
            },
            index=[0],
        )
        return float(np.linalg.norm(demand.calculate_forces_vectorized(particles)[0]))

    capillary = pull(0.004)
    assert capillary > 0
    np.testing.assert_allclose(pull(0.016), capillary * 0.25, rtol=1e-9)
    # At or below the reference caliber the pull is unattenuated
    np.testing.assert_allclose(pull(0.002), capillary, rtol=1e-9)


def test_perfusion_demand_improves_tissue_coverage():
    # Deterministic seeds make this an exact comparison, not a statistical one
    coverage_without = run_simulation(with_perfusion=False)
    coverage_with = run_simulation(with_perfusion=True)
    assert coverage_with > coverage_without, (
        f"perfusion demand should improve coverage: "
        f"{coverage_with:.3f} vs {coverage_without:.3f}"
    )

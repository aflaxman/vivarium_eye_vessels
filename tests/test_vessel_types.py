"""Tests for paired arterial and venous trees (roadmap idea 3)."""

import copy

import numpy as np
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import (
    EllipsoidContainment,
    FrozenRepulsion,
    PerfusionDemand,
)
from vivarium_eye_vessels.components.particles import (
    VESSEL_TYPE_ARTERY,
    VESSEL_TYPE_NONE,
    VESSEL_TYPE_VEIN,
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
        "artery_caliber_ratio": 0.67,
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
        "cross_type_factor": 0.25,
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


def make_simulation() -> InteractiveContext:
    return InteractiveContext(
        components=[
            Particle3D(),
            PathFreezer(),
            PathExtinction(),
            PathSplitter(),
            EllipsoidContainment(),
            FrozenRepulsion(),
            PerfusionDemand(),
        ],
        configuration=copy.deepcopy(CONFIGURATION),
    )


def test_roots_alternate_types_with_avr_calibers():
    sim = make_simulation()
    pop = sim.get_population(["vessel_type", "path_id", "radius"])
    roots = pop[pop.path_id >= 0].sort_index()

    assert list(roots.vessel_type) == [
        VESSEL_TYPE_ARTERY,
        VESSEL_TYPE_VEIN,
        VESSEL_TYPE_ARTERY,
        VESSEL_TYPE_VEIN,
    ]
    artery_radius = roots[roots.vessel_type == VESSEL_TYPE_ARTERY].radius.iloc[0]
    vein_radius = roots[roots.vessel_type == VESSEL_TYPE_VEIN].radius.iloc[0]
    np.testing.assert_allclose(artery_radius / vein_radius, 0.67, atol=1e-12)

    off_path = pop[pop.path_id < 0]
    assert (off_path.vessel_type == VESSEL_TYPE_NONE).all()


def test_types_are_inherited_down_both_trees():
    sim = make_simulation()
    simulation.run_steps(sim, 150)
    pop = sim.get_population(["vessel_type", "path_id", "parent_id", "frozen"])

    on_path = pop[pop.path_id >= 0]
    assert set(on_path.vessel_type.unique()) == {VESSEL_TYPE_ARTERY, VESSEL_TYPE_VEIN}
    assert (on_path.vessel_type == VESSEL_TYPE_ARTERY).sum() > 10
    assert (on_path.vessel_type == VESSEL_TYPE_VEIN).sum() > 10

    # Every child carries its parent's type
    children = on_path[on_path.parent_id >= 0]
    children = children[children.parent_id.isin(pop.index)]
    parent_types = pop.loc[children.parent_id, "vessel_type"].values
    assert (children.vessel_type.values == parent_types).all()


def test_weak_cross_type_repulsion_lets_trees_interdigitate():
    """With weaker artery-vein repulsion, the trees end up closer together.

    Deterministic seeds make this an exact comparison, not a statistical one.
    """

    def median_cross_tree_distance(cross_type_factor: float) -> float:
        from scipy.spatial import cKDTree

        config = copy.deepcopy(CONFIGURATION)
        config["frozen_repulsion"]["cross_type_factor"] = cross_type_factor
        sim = InteractiveContext(
            components=[
                Particle3D(),
                PathFreezer(),
                PathExtinction(),
                PathSplitter(),
                EllipsoidContainment(),
                FrozenRepulsion(),
                PerfusionDemand(),
            ],
            configuration=config,
        )
        simulation.run_steps(sim, 250)
        pop = sim.get_population(["x", "y", "z", "vessel_type", "frozen"])
        arteries = pop[pop.frozen & (pop.vessel_type == VESSEL_TYPE_ARTERY)]
        veins = pop[pop.frozen & (pop.vessel_type == VESSEL_TYPE_VEIN)]
        assert len(arteries) > 10 and len(veins) > 10
        distances, _ = cKDTree(veins[["x", "y", "z"]].to_numpy()).query(
            arteries[["x", "y", "z"]].to_numpy(), k=1
        )
        return float(np.median(distances))

    assert median_cross_tree_distance(0.25) < median_cross_tree_distance(1.0)


def test_typed_perfusion_demand_keeps_trees_balanced():
    """Per-type demand means neither tree can win territory for both."""
    sim = make_simulation()
    simulation.run_steps(sim, 250)
    pop = sim.get_population(["vessel_type", "frozen", "path_id"])
    frozen = pop[pop.frozen & (pop.path_id >= 0)]
    n_artery = int((frozen.vessel_type == VESSEL_TYPE_ARTERY).sum())
    n_vein = int((frozen.vessel_type == VESSEL_TYPE_VEIN).sum())
    assert n_artery > 50 and n_vein > 50
    balance = min(n_artery, n_vein) / max(n_artery, n_vein)
    assert balance > 0.4, f"trees out of balance: {n_artery} arteries vs {n_vein} veins"

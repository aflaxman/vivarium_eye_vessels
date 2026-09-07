"""Tests for the capillary bed: fine sites, sprouting from vessel walls, starving, and the hooks."""

import copy

import numpy as np
import pandas as pd
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import (
    CapillaryBed,
    CylinderExclusion,
    EllipsoidContainment,
    FrozenRepulsion,
    PerfusionDemand,
    capillary_sites,
)
from vivarium_eye_vessels.components.particles import (
    Particle3D,
    PathExtinction,
    PathFreezer,
    PathSplitter,
)

CONFIGURATION = {
    "randomness": {"random_seed": 3},
    "time": {
        "start": {"year": 2025, "month": 1, "day": 1},
        "end": {"year": 2025, "month": 3, "day": 1},
        "step_size": 0.05,
    },
    "population": {"population_size": 400},
    "particles": {
        "initial_circle": {"center": [0.5, 0.0, 0.0], "radius": 0.05, "n_vessels": 4},
        "root_radius": 0.017,
        "terminal_velocity": 0.18,
    },
    "path_splitter": {"split_interval": 5, "min_radius": 0.002, "capillary_radius": 0.00095},
    "ellipsoid_containment": {"a": 1.0, "b": 1.0, "c": 0.1, "force_type": "hookean"},
    "perfusion_demand": {"site_spacing": 0.1, "perfusion_radius": 0.15},
    "cylinder_exclusion": {
        "radius": 0.05,
        "center": [-0.2, 0.0, 0.0],
        "direction": [0.0, 0.0, 1.0],
        "force_type": "hookean",
        "spring_constant": 30,
    },
    "frozen_repulsion": {
        "force_type": "hookean",
        "spring_constant": 1.5,
        "interaction_radius": 0.12,
        "capillary_radius": 0.0015,
        "capillary_interaction_radius": 0.02,
    },
    "capillary_bed": {
        "enabled": True,
        "region_radius": 0.0,
        "site_spacing": 0.05,
        "perfusion_radius": 0.03,
        "capillary_radius": 0.0009,
        "sprout_interval": 1,
        "sprout_range": 0.1,
        "max_sprouts": 10,
        "starve_radius": 0.04,
    },
}


def build(steps: int = 12):
    bed, repulsion, splitter = CapillaryBed(), FrozenRepulsion(), PathSplitter()
    sim = InteractiveContext(
        components=[
            Particle3D(),
            PathFreezer(),
            PathExtinction(),
            splitter,
            EllipsoidContainment(),
            CylinderExclusion(),
            repulsion,
            PerfusionDemand(),
            bed,
        ],
        configuration=copy.deepcopy(CONFIGURATION),
    )
    for _ in range(steps):
        sim.step()
    return sim, bed, repulsion, splitter


def population(sim, bed):
    return sim.get_population(bed.required_attributes)


def test_capillary_sites_cover_the_region_in_every_layer_and_skip_the_fovea():
    positions, layers = capillary_sites(
        (1.0, 1.0), 0.05, (-0.2, 0.0), 0.4, 0.05, [0.04, 0.0, -0.04]
    )
    assert set(np.unique(layers)) == {0, 1, 2}
    per_layer = np.bincount(layers)
    assert per_layer[0] == per_layer[1] == per_layer[2] > 100
    from_fovea = np.hypot(positions[:, 0] + 0.2, positions[:, 1])
    assert from_fovea.max() <= 0.4 and from_fovea.min() > 0.05
    np.testing.assert_allclose(np.unique(positions[:, 2]), [-0.04, 0.0, 0.04])


def test_capillary_tips_never_split():
    sim, bed, _, splitter = build(steps=1)
    tips = pd.DataFrame({"radius": [0.0009, 0.002, 0.005], "vessel_type": [1, 1, 1]})
    probabilities = splitter.split_probabilities(tips)
    assert probabilities.iloc[0] == 0.0 and probabilities.iloc[1] > 0


def test_sprouts_are_capillary_tips_aimed_from_a_wall_at_hypoxic_tissue():
    sim, bed, _, splitter = build()
    pop = population(sim, bed)
    before = set(pop.index[~pop.frozen & (pop.path_id >= 0)])
    frozen = pop[pop.frozen & (pop.path_id >= 0)]
    assert not frozen.empty
    sites, layers = bed.hypoxic_sites()
    assert len(sites) > 0
    next_id = splitter.next_path_id
    bed.sprout(pop, sites, layers)
    after = population(sim, bed)
    sprouts = after[~after.frozen & (after.path_id >= 0) & ~after.index.isin(list(before))]
    assert 0 < len(sprouts) <= 10
    np.testing.assert_allclose(sprouts.radius, 0.0009)
    assert sprouts.parent_id.isin(frozen.index).all()
    assert sprouts.path_id.min() >= next_id and sprouts.path_id.is_unique
    # Each sprout heads straight from its wall at one of the hypoxic sites in range
    origins = pop.loc[sprouts.parent_id, ["x", "y", "z"]].to_numpy()
    headings = sprouts[["vx", "vy", "vz"]].to_numpy()
    headings /= np.linalg.norm(headings, axis=1)[:, np.newaxis]
    for origin, heading in zip(origins, headings):
        offsets = sites - origin
        distances = np.linalg.norm(offsets, axis=1)
        in_range = distances <= 0.1
        cosines = offsets[in_range] @ heading / distances[in_range]
        assert cosines.max() > 0.99


def test_capillary_tips_starve_where_no_tissue_calls_them():
    sim, bed, _, _ = build()
    pop = population(sim, bed)
    sites, layers = bed.hypoxic_sites()
    bed.sprout(pop, sites, layers)
    pop = population(sim, bed)
    capillaries = pop[~pop.frozen & (pop.path_id >= 0) & (pop.radius <= 0.0009)]
    assert not capillaries.empty
    # With no hypoxic site anywhere, every capillary tip freezes
    bed.starve(pop, np.empty((0, 3)), np.empty(0, dtype=int))
    pop = population(sim, bed)
    assert pop.loc[capillaries.index, "frozen"].all()
    assert (pop.loc[capillaries.index, "path_id"] == -1).all()


def test_capillary_tips_feel_only_short_range_repulsion():
    sim, bed, repulsion, _ = build()
    pop = population(sim, bed)
    wall = pop[pop.frozen & (pop.path_id >= 0)].iloc[0]
    tips = pd.DataFrame(
        {
            "x": wall.x + 0.05,
            "y": wall.y,
            "z": wall.z,
            "frozen": False,
            "path_id": [9001, 9002],
            "parent_id": -1,
            "freeze_time": pd.NaT,
            "vessel_type": 1,
            "radius": [0.001, 0.005],
        }
    )
    forces = repulsion.calculate_forces_vectorized(tips)
    assert np.linalg.norm(forces[0]) == 0.0  # a capillary 225 um away is out of its reach
    assert np.linalg.norm(forces[1]) > 0.0  # an arteriole tip is still repelled


def test_capillaries_do_not_perfuse_the_coarse_lattice():
    sim, bed, _, _ = build()
    demand = [c for c in sim.list_components().values() if isinstance(c, PerfusionDemand)][0]
    demand.min_radius = 0.00095
    frozen = bed.freezer.frozen_particles()
    site = frozen[["x", "y", "z"]].to_numpy()[:1]
    with_all = demand.vessel_distances(site)
    demand.min_radius = 0.0
    assert with_all is not None and np.isclose(with_all[0], 0.0)
    # Turn the nearest frozen vessel into a capillary: the site is no longer served by it
    bed.particles.update_particles(pd.DataFrame({"radius": 0.0009}, index=frozen.index[:1]))
    bed.freezer.update_tree(sim.get_population(bed.freezer.required_attributes))
    demand.min_radius = 0.00095
    assert demand.vessel_distances(site)[0] > 0.0


def test_wide_tips_are_not_fenced_out_by_frozen_capillaries():
    sim, bed, repulsion, _ = build()
    bed.freezer.update_tree(sim.get_population(bed.freezer.required_attributes))
    frozen = bed.freezer.frozen_particles()
    wall = frozen.iloc[0]
    tip = pd.DataFrame(
        {
            "x": [wall.x + 0.03],
            "y": [wall.y],
            "z": [wall.z],
            "frozen": False,
            "path_id": [9001],
            "parent_id": -1,
            "freeze_time": pd.NaT,
            "vessel_type": 1,
            "radius": [0.005],
        }
    )
    before = np.linalg.norm(repulsion.calculate_forces_vectorized(tip)[0])
    assert before > 0
    # Turn every frozen vessel within reach into a capillary: the arteriole tip feels nothing
    distance = np.sqrt(
        (frozen.x - tip.x[0]) ** 2 + (frozen.y - tip.y[0]) ** 2 + (frozen.z - tip.z[0]) ** 2
    )
    near = frozen[distance < 0.12].index
    bed.particles.update_particles(pd.DataFrame({"radius": 0.0009}, index=near))
    bed.freezer.update_tree(sim.get_population(bed.freezer.required_attributes))
    assert np.linalg.norm(repulsion.calculate_forces_vectorized(tip)[0]) == 0.0


def test_precapillary_tips_do_not_end_in_the_capillary_bed():
    from vivarium_eye_vessels.components.particles import anastomosis_targets

    tips = pd.DataFrame(
        {
            "x": [0.0, 0.0],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0],
            "vessel_type": [1, 1],
            "radius": [0.003, 0.0009],
        },
        index=[10, 11],
    )
    frozen = pd.DataFrame(
        {
            "x": [0.01, 0.02],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0],
            "vessel_type": [2, 2],
            "radius": [0.0009, 0.003],
            "layer_id": [0, 0],
        },
        index=[20, 21],
    )
    neighbors = [[0, 1], [0, 1]]
    plain = anastomosis_targets(tips, frozen, neighbors, max_target_radius=0.004)
    assert plain[10] == 20 and plain[11] == 20  # nearest opposite-tree segment, a capillary
    guarded = anastomosis_targets(
        tips, frozen, neighbors, max_target_radius=0.004, capillary_radius=0.00095
    )
    assert guarded[10] == 21  # the precapillary tip skips the capillary for the venule
    assert guarded[11] == 20  # the capillary tip still closes its loop on a capillary


def test_dead_end_sprouts_regress_and_closed_loops_stay():
    sim, bed, _, _ = build()
    pop = population(sim, bed)
    sites, layers = bed.hypoxic_sites()
    bed.sprout(pop, sites, layers)
    pop = population(sim, bed)
    tips = pop[~pop.frozen & (pop.path_id >= 0) & (pop.radius <= 0.0009)]
    assert len(tips) >= 2
    # Stop every sprout where it stands (dead ends), then mark one as a closed loop
    bed.starve(pop, np.empty((0, 3)), np.empty(0, dtype=int))
    keeper = tips.index[0]
    bed.particles.update_particles(
        pd.DataFrame({"anastomosis_id": [int(tips.parent_id.iloc[1])]}, index=[keeper])
    )
    # Too young to regress
    pop = population(sim, bed)
    bed.regress(pop)
    pop = population(sim, bed)
    assert pop.loc[tips.index, "frozen"].all()
    # Old enough: the dead ends are recycled, the closed loop and the walls survive
    walls = set(tips.parent_id)
    bed.particles.update_particles(
        pd.DataFrame({"freeze_time": bed.clock() - pd.Timedelta(days=10)}, index=tips.index)
    )
    pop = population(sim, bed)
    bed.regress(pop)
    pop = population(sim, bed)
    others = tips.index.drop(keeper)
    assert (pop.loc[others, "path_id"] == -1).all() and (~pop.loc[others, "frozen"]).all()
    assert (pop.loc[others, "radius"] == 0).all()
    assert pop.loc[keeper, "frozen"] and pop.loc[keeper, "radius"] == 0.0009
    # The walls the sprouts left (arterioles, or older capillaries still within
    # their grace) stand where they were
    assert pop.loc[list(walls), "frozen"].all()
    arteriole_walls = [w for w in walls if pop.loc[w, "radius"] > 0.001]
    assert arteriole_walls and (pop.loc[arteriole_walls, "path_id"] >= 0).all()


def test_a_capillary_tip_inside_the_fovea_is_withdrawn_not_frozen_there():
    sim, bed, _, _ = build()
    pop = population(sim, bed)
    sites, layers = bed.hypoxic_sites()
    bed.sprout(pop, sites, layers)
    pop = population(sim, bed)
    tip = pop[~pop.frozen & (pop.path_id >= 0) & (pop.radius <= 0.0009)].index[0]
    # Move the tip into the avascular zone (fovea at (-0.2, 0), radius 0.05)
    bed.particles.update_particles(pd.DataFrame({"x": [-0.2], "y": [0.01]}, index=[tip]))
    bed.starve(population(sim, bed), sites, layers)
    pop = population(sim, bed)
    assert (
        not pop.loc[tip, "frozen"]
        and pop.loc[tip, "path_id"] == -1
        and pop.loc[tip, "radius"] == 0
    )

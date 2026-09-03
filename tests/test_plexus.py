"""Tests for stratified vascular plexuses (roadmap idea 6)."""

import copy

import numpy as np
from vivarium import InteractiveContext

from vivarium_eye_vessels.components.boundaries import (
    EllipsoidContainment,
    FrozenRepulsion,
    PerfusionDemand,
)
from vivarium_eye_vessels.components.particles import (
    Particle3D,
    PathExtinction,
    PathFreezer,
    PathSplitter,
)
from vivarium_eye_vessels.components.plexus import PlexusLayers
from vivarium_eye_vessels.vnv import simulation

LAYER_Z = [0.04, 0.0, -0.04]

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
        # Roots start on the superficial plane, so the stratification test
        # measures how well vessels hold their layer, not how fast the
        # spring pulls an off-plane start onto it
        "initial_circle": {"center": [1.0, 0.0, LAYER_Z[0]], "radius": 0.05, "n_vessels": 4},
        "root_radius": 0.02,
        "artery_caliber_ratio": 0.67,
    },
    "path_freezer": {"freeze_interval": 3, "radius_taper": 0.99},
    "path_splitter": {
        "split_interval": 15,
        "split_angle": 60,
        "split_probability": 0.9,
        "max_depth": 4,
        "murray_exponent": 3.0,
        "flow_asymmetry": 0.15,
        "min_radius": 0.002,
        "caliber_cadence_exponent": 0.0,
    },
    "path_extinction": {"force_threshold": 1.2},
    "plexus_layers": {
        "layer_z": LAYER_Z,
        "spring_constant": 6.0,
        "damping": 5.0,
        "max_force": 0.3,
        "dive_radius": 0.01,
        "dive_probability": 0.05,
    },
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
            PlexusLayers(),
            EllipsoidContainment(),
            FrozenRepulsion(),
            PerfusionDemand(),
        ],
        configuration=copy.deepcopy(CONFIGURATION),
    )


def test_roots_start_superficial_and_layers_get_populated():
    sim = make_simulation()
    pop = sim.get_population(["layer_id", "path_id"])
    roots = pop[pop.path_id >= 0]
    assert (roots.layer_id == 0).all()
    assert (pop[pop.path_id < 0].layer_id == -1).all()

    simulation.run_steps(sim, 250)
    pop = sim.get_population(["layer_id", "path_id", "frozen", "radius"])
    vessels = pop[(pop.layer_id >= 0) & pop.frozen]
    populated = set(vessels.layer_id.unique())
    assert 0 in populated
    assert populated - {0}, "no tips ever dove below the superficial plexus"


def test_only_capillary_tips_dive():
    sim = make_simulation()
    simulation.run_steps(sim, 250)
    pop = sim.get_population(["layer_id", "radius", "frozen", "path_id", "parent_id"])
    vessels = pop[(pop.layer_id >= 1) & (pop.radius > 0)]
    assert len(vessels) > 0
    dive_radius = CONFIGURATION["plexus_layers"]["dive_radius"]
    taper = CONFIGURATION["path_freezer"]["radius_taper"]
    # Deep-layer vessels descend from capillary-caliber divers; taper only
    # shrinks calibers, so nothing wider than the dive threshold can be there
    assert vessels.radius.max() <= dive_radius / taper + 1e-12


def test_vessels_stratify_onto_their_planes():
    sim = make_simulation()
    simulation.run_steps(sim, 250)
    pop = sim.get_population(["z", "layer_id", "frozen", "path_id"])
    vessels = pop[(pop.layer_id >= 0) & pop.frozen]

    planes = np.asarray(LAYER_Z)
    z_errors = (vessels.z - planes[vessels.layer_id.to_numpy(int)]).abs()
    layer_gap = abs(LAYER_Z[0] - LAYER_Z[1])
    # Vessels sit closer to their own plane than to the next plexus over
    assert float(z_errors.median()) < layer_gap / 2


def test_diving_vessels_connect_the_layers():
    sim = make_simulation()
    simulation.run_steps(sim, 250)
    pop = sim.get_population(["layer_id", "parent_id", "path_id", "frozen"])
    on_path = pop[pop.layer_id >= 0]
    children = on_path[on_path.parent_id.isin(on_path.index)]
    parent_layers = on_path.layer_id.loc[children.parent_id].to_numpy()
    dives = children.layer_id.to_numpy() != parent_layers
    assert dives.sum() > 0, "no diving vessels formed"
    # Diving only goes downward (a tip may dive more than once between
    # freezes, so a recorded parent-child hop can span two layers)
    deltas = (children.layer_id.to_numpy() - parent_layers)[dives]
    assert (deltas >= 1).all()
    assert deltas.max() <= len(LAYER_Z) - 1

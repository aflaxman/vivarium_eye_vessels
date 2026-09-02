from typing import List

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import norm
from vivarium import Component
from vivarium.framework.engine import Builder
from vivarium.framework.event import Event
from vivarium.framework.population import SimulantData

PARTICLE_COLUMNS = [
    # location
    "x",
    "y",
    "z",
    # velocity
    "vx",
    "vy",
    "vz",
    # "freeze" information used to form eye vessels
    "frozen",
    "freeze_time",
    "unfreeze_time",
    "depth",
    # addl information relevant to eye vessel structure
    "parent_id",  # tree structure
    "path_id",  # used to hack PathExtinction dynamics so that splits don't go extinct immediately
    "radius",  # vessel caliber, assigned by Murray's law at bifurcations
    "vessel_type",  # artery/vein identity, inherited down each tree
    "anastomosis_id",  # index of the opposite-tree particle this tip fused with
    "layer_id",  # which vascular plexus the vessel belongs to (0 = superficial)
    # Ornstein-Uhlenbeck steering state: the autocorrelated random component
    # of each tip's acceleration (see particles.noise_persistence_time)
    "wx",
    "wy",
    "wz",
]

VESSEL_TYPE_NONE = 0
VESSEL_TYPE_ARTERY = 1
VESSEL_TYPE_VEIN = 2


def murray_daughter_radii(parent_radius: float, flow_fraction, exponent: float) -> tuple:
    """Daughter radii for a bifurcation obeying Murray's law.

    Flow splits into ``flow_fraction`` (minor daughter) and ``1 - flow_fraction``
    (major daughter); with flow proportional to radius**exponent, the radii
    satisfy r_minor**k + r_major**k == r_parent**k.
    """
    minor = parent_radius * flow_fraction ** (1.0 / exponent)
    major = parent_radius * (1.0 - flow_fraction) ** (1.0 / exponent)
    return major, minor


def murray_bifurcation_angles(parent_radius, major_radius, minor_radius) -> tuple:
    """Optimal daughter angles (radians, relative to the parent axis).

    Murray's minimum-work principle (Murray 1926) gives, for parent radius r0
    and daughters r1, r2:

        cos(theta_1) = (r0**4 + r1**4 - r2**4) / (2 r0**2 r1**2)
        cos(theta_2) = (r0**4 + r2**4 - r1**4) / (2 r0**2 r2**2)

    so the large daughter deviates little from the parent axis and the small
    daughter comes off at a wider angle.
    """
    r0, r1, r2 = parent_radius, major_radius, minor_radius
    cos_1 = (r0**4 + r1**4 - r2**4) / (2 * r0**2 * r1**2)
    cos_2 = (r0**4 + r2**4 - r1**4) / (2 * r0**2 * r2**2)
    return np.arccos(np.clip(cos_1, -1.0, 1.0)), np.arccos(np.clip(cos_2, -1.0, 1.0))


class Particle3D(Component):
    """Base component for managing 3D particle positions, velocities, and forces.

    Under vivarium 4, only the component that creates a column may write to it,
    so this component owns all particle state and exposes :meth:`update_particles`
    for sibling components (PathFreezer, PathSplitter, etc.) to route their
    state changes through.
    """

    CONFIGURATION_DEFAULTS = {
        "particles": {
            "overall_max_velocity_change": 0.1,
            "initial_velocity_range": (-0.05, 0.05),
            "terminal_velocity": 0.2,  # Maximum allowed velocity magnitude
            "initial_circle": {"center": [1.5, 0.0, 0.5], "radius": 0.1, "n_vessels": 5},
            "root_radius": 0.02,  # Caliber of the root (vein) vessels at the seed circle
            # Arteries are narrower than veins; the clinical artery:vein ratio is ~2:3
            "artery_caliber_ratio": 0.67,
            # Persistence time (days) of the random steering. 0 keeps the
            # legacy white-noise kicks; above 0 the kicks become an
            # Ornstein-Uhlenbeck process with this correlation time, whose
            # coherent curvature raises tortuosity relative to white noise --
            # the disease dial (see realism roadmap idea 7)
            "noise_persistence_time": 0.0,
            # Caliber-dependent steering stiffness: a tip wider than the
            # reference caliber has its random steering attenuated by
            # (reference / radius) ** exponent, so arcade-caliber tips hold
            # their heading while capillary tips wander freely. Exponent 0
            # keeps the legacy caliber-blind steering
            "noise_caliber_reference": 0.004,
            "noise_caliber_exponent": 0.0,
        }
    }

    def setup(self, builder: Builder) -> None:
        self.config = builder.configuration.particles
        self.step_size = builder.configuration.time.step_size
        self.overall_max_velocity_change = self.config.overall_max_velocity_change
        self.initial_velocity_range = self.config.initial_velocity_range
        self.terminal_velocity = self.config.terminal_velocity
        self.noise_persistence_time = float(self.config.noise_persistence_time)
        self.noise_caliber_reference = float(self.config.noise_caliber_reference)
        self.noise_caliber_exponent = float(self.config.noise_caliber_exponent)

        self.clock = builder.time.clock()

        # Register force pipelines
        self.register_force_pipelines(builder)

        self.max_velocity_change = builder.value.register_value_producer(
            "particle.max_velocity_change",
            source=lambda index: pd.Series(self.overall_max_velocity_change, index=index),
        )

        self.randomness = builder.randomness.get_stream("particle.particles_3d")
        builder.population.register_initializer(
            initializer=self.on_initialize_simulants,
            columns=PARTICLE_COLUMNS,
            required_resources=[self.randomness],
        )
        self.setup_scale(builder)

    def setup_scale(self, builder: Builder) -> None:
        has_ellipsoid = "ellipsoid_containment" in builder.components.list_components()

        if has_ellipsoid:
            # Get ellipsoid parameters
            config = builder.configuration.ellipsoid_containment
            a = float(config.a)
            b = float(config.b)
            c = float(config.c)
            self.scale = np.array([a, b, c])
        else:
            self.scale = np.ones(3)

    def register_force_pipelines(self, builder: Builder) -> None:
        """Register pipelines for force components and total magnitude."""
        # Register individual force component pipelines
        self.force_x = builder.value.register_value_producer(
            "particle.force.x", source=lambda index: pd.Series(0.0, index=index)
        )
        self.force_y = builder.value.register_value_producer(
            "particle.force.y", source=lambda index: pd.Series(0.0, index=index)
        )
        self.force_z = builder.value.register_value_producer(
            "particle.force.z", source=lambda index: pd.Series(0.0, index=index)
        )

        # Register total force magnitude pipeline
        self.force_magnitude = builder.value.register_value_producer(
            "particle.force.magnitude",
            source=self.get_force_magnitude,
            required_resources=[self.force_x, self.force_y, self.force_z],
        )

    def get_force_magnitude(self, index: pd.Index) -> pd.Series:
        """Calculate total force magnitude from components."""
        fx = self.force_x(index)
        fy = self.force_y(index)
        fz = self.force_z(index)
        return np.sqrt(fx**2 + fy**2 + fz**2)

    def get_particles(self, index: pd.Index) -> pd.DataFrame:
        """Get the full particle state table for the given index."""
        return self.population_view.get(index, PARTICLE_COLUMNS)

    def update_particles(self, updates: pd.DataFrame) -> None:
        """Write particle state updates on behalf of this or a sibling component.

        Only this component may write the particle columns, so components that
        change particle state (freezing, splitting, extinction) construct a
        DataFrame of new values indexed by simulant and pass it here.
        """
        # Coerce datetime columns to the state table's resolution (pandas 3
        # infers datetime64[us] for new frames, but the table uses ns).
        datetime_columns = [
            col for col in ("freeze_time", "unfreeze_time") if col in updates.columns
        ]
        if datetime_columns:
            updates = updates.copy()
            for col in datetime_columns:
                updates[col] = updates[col].astype("datetime64[ns]")
        columns = list(updates.columns)
        self.population_view.update(columns, lambda current: updates)

    def on_initialize_simulants(self, simulant_data: SimulantData) -> None:
        """Initialize particles with positions, velocities, and path tracking information."""
        pop = pd.DataFrame(index=simulant_data.index)

        # Generate 3D normal points using ppf (inverse CDF)
        points = np.column_stack(
            [
                norm.ppf(self.randomness.get_draw(pop.index, additional_key=f"xyz_{i}"))
                for i in range(3)
            ]
        )

        # Normalize and scale by random radius
        points /= np.linalg.norm(points, axis=1)[:, np.newaxis]
        radii = np.array(self.randomness.get_draw(pop.index, additional_key="radius")) ** (
            1 / 3
        )
        points *= radii[:, np.newaxis]
        pop[["x", "y", "z"]] = points * self.scale

        # Generate random initial velocities
        v_range = self.initial_velocity_range
        for v in ["vx", "vy", "vz"]:
            pop[v] = (
                self.randomness.get_draw(pop.index, additional_key=v)
                * (v_range[1] - v_range[0])
                + v_range[0]
            )
        pop[["vx", "vy", "vz"]] *= self.scale

        # Initialize tree-structure-related columns
        pop["frozen"] = False
        pop["freeze_time"] = pd.NaT
        pop["unfreeze_time"] = pd.NaT
        pop["depth"] = -1
        pop["parent_id"] = -1
        pop["path_id"] = -1
        pop["radius"] = 0.0
        pop["vessel_type"] = VESSEL_TYPE_NONE
        pop["anastomosis_id"] = -1
        pop["layer_id"] = -1
        pop["wx"] = 0.0
        pop["wy"] = 0.0
        pop["wz"] = 0.0

        self.initialize_circle_positions(pop)

        self.population_view.initialize(pop)

    def initialize_circle_positions(self, pop: pd.DataFrame) -> None:
        # Initialize active vessel in circle position
        config = self.config.initial_circle
        center = config.center
        radius = config.radius
        n_vessels = config.n_vessels

        for i in range(n_vessels):
            if i in pop.index:
                angle = 2 * np.pi * i / n_vessels
                pop.loc[i, ["x", "y", "z"]] = [
                    center[0] + radius * np.cos(angle),
                    center[1] + radius * np.sin(angle),
                    center[2],
                ]
                pop.loc[i, ["vx", "vy", "vz"]] = [
                    radius * np.cos(angle),
                    radius * np.sin(angle),
                    0.0,
                ]
                pop.loc[i, "path_id"] = i
                pop.loc[i, ["depth"]] = 0
                pop.loc[i, "layer_id"] = 0  # roots live in the superficial plexus
                # Alternate artery/vein arcades around the disc; arteries are narrower
                if i % 2 == 0:
                    pop.loc[i, "vessel_type"] = VESSEL_TYPE_ARTERY
                    pop.loc[i, "radius"] = (
                        self.config.root_radius * self.config.artery_caliber_ratio
                    )
                else:
                    pop.loc[i, "vessel_type"] = VESSEL_TYPE_VEIN
                    pop.loc[i, "radius"] = self.config.root_radius

    def on_time_step(self, event: Event) -> None:
        """Update positions and velocities of non-frozen particles and track blocking forces."""
        pop = self.get_particles(event.index)
        active_particles = pop[~pop.frozen]

        if not active_particles.empty:
            self.update_positions(active_particles)

    def update_positions(self, particles: pd.DataFrame) -> None:
        """Update positions and velocities based on forces and random steering."""
        columns = ["x", "y", "z", "vx", "vy", "vz"]
        if self.noise_persistence_time > 0:
            columns += ["wx", "wy", "wz"]
        updates = particles[columns].copy()

        # Update positions based on current velocities
        for pos, vel in [("x", "vx"), ("y", "vy"), ("z", "vz")]:
            updates[pos] = updates[pos] + self.step_size * updates[vel]

        # Get max velocity change from pipeline
        max_velocity_change = self.max_velocity_change(updates.index)

        # Get current forces from pipelines
        fx = self.force_x(updates.index)
        fy = self.force_y(updates.index)
        fz = self.force_z(updates.index)

        # Wide tips steer stiffly: attenuate the random kick with caliber
        if self.noise_caliber_exponent > 0:
            radius = particles.radius.to_numpy()
            attenuation = np.where(
                radius > self.noise_caliber_reference,
                (self.noise_caliber_reference / np.maximum(radius, 1e-12))
                ** self.noise_caliber_exponent,
                1.0,
            )
        else:
            attenuation = 1.0

        # Update velocities with random steering and forces
        for i, (v, w, f) in enumerate(
            zip(["vx", "vy", "vz"], ["wx", "wy", "wz"], [fx, fy, fz])
        ):
            if self.noise_persistence_time > 0:
                # Ornstein-Uhlenbeck steering: an AR(1) with correlation time
                # noise_persistence_time and the same stationary spread as the
                # legacy uniform kick (sd = max_velocity_change / sqrt(3)), so
                # theta -> 1 degenerates exactly to white noise
                theta = min(self.step_size / self.noise_persistence_time, 1.0)
                stationary_sd = (
                    max_velocity_change / np.sqrt(3.0) * self.scale[i] * attenuation
                )
                draws = self.randomness.get_draw(updates.index, additional_key=f"d{v}")
                shocks = norm.ppf(np.clip(draws, 1e-12, 1 - 1e-12))
                updates[w] = (1 - theta) * updates[w] + stationary_sd * np.sqrt(
                    theta * (2 - theta)
                ) * shocks
                dv = updates[w]
            else:
                # Legacy white-noise kick, uniform on +/- max_velocity_change
                dv = (
                    (self.randomness.get_draw(updates.index, additional_key=f"d{v}") - 0.5)
                    * 2
                    * max_velocity_change
                    * self.scale[i]
                    * attenuation
                )

            # Add force contribution to velocity
            updates[v] += (dv + f) * self.step_size

        # Apply terminal velocity constraint
        velocity_vectors = updates[["vx", "vy", "vz"]].to_numpy() / self.scale
        velocities_magnitude = np.linalg.norm(velocity_vectors, axis=1)
        over_limit = velocities_magnitude > self.terminal_velocity

        if np.any(over_limit):
            # Scale down velocity components to satisfy terminal velocity
            scale_factors = self.terminal_velocity / velocities_magnitude[over_limit]
            updates.loc[over_limit, ["vx", "vy", "vz"]] *= scale_factors[:, np.newaxis]

        self.update_particles(updates)


class PathFreezer(Component):
    """Component for freezing particle paths and creating continuations.
    Also mantains KDTree of frozen particles for efficient querying."""

    CONFIGURATION_DEFAULTS = {
        "path_freezer": {
            "freeze_interval": 10,
            "radius_taper": 1.0,  # caliber multiplier per frozen segment along a path
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return [
            "x",
            "y",
            "z",
            "vx",
            "vy",
            "vz",
            "frozen",
            "freeze_time",
            "depth",
            "parent_id",
            "path_id",
            "radius",
            "vessel_type",
            "layer_id",
        ]

    def setup(self, builder: Builder) -> None:
        self.config = builder.configuration.path_freezer
        self.particles_to_add = builder.configuration.population.population_size
        self.step_count = 0
        self.clock = builder.time.clock()

        self._current_tree = None
        self._current_frozen = None
        self.simulant_creator = builder.population.get_simulant_creator()
        self.particles = builder.components.get_components_by_type(Particle3D)[0]

    def add_particles(self):
        self.simulant_creator(self.particles_to_add)

    def on_time_step(self, event: Event) -> None:
        self.step_count += 1
        if self.step_count % self.config.freeze_interval == 0:
            pop = self.population_view.get(event.index, self.required_attributes)
            self.freeze_particles(pop)
            self.update_tree(pop)

    def update_tree(self, pop):
        self._current_frozen = pop[pop.frozen]
        if len(self._current_frozen) < 2:
            self._current_tree = None
        else:
            self._current_tree = cKDTree(self._current_frozen[["x", "y", "z"]].values)

    def frozen_particles(self):
        """The frozen-particle table backing the KDTree, or None before any freeze."""
        return self._current_frozen

    def get_neighbor_pairs(self, radius: float):
        """Get all pairs of frozen particles within radius using efficient pair query."""
        if self._current_tree is None:
            return None

        return self._current_tree.query_pairs(radius)

    def query_radius(self, pop, radius: float):
        """Get neighbor indices for each particle within radius."""
        if self._current_tree is None:
            return None

        if isinstance(pop, pd.DataFrame):
            positions = pop[["x", "y", "z"]].values
        else:
            positions = pop

        return self._current_tree.query_ball_point(positions, radius)

    def get_population(self, indices: List[int]) -> pd.DataFrame:
        """Get frozen particles by their positional indices in the KDTree."""
        return self._current_frozen.iloc[list(indices)]

    def freeze_particles(self, pop: pd.DataFrame) -> None:
        """Create frozen path points and continue paths with new particles."""
        active = pop[~pop.frozen & (pop.path_id >= 0)]
        if active.empty:
            return

        available = pop[~pop.frozen & (pop.path_id < 0)]
        if len(available) >= len(active):
            to_freeze = available.iloc[: len(active)]

            continuations = pd.DataFrame(
                {
                    "x": active.x.values,
                    "y": active.y.values,
                    "z": active.z.values,
                    "vx": active.vx.values,
                    "vy": active.vy.values,
                    "vz": active.vz.values,
                    "path_id": active.path_id.values,
                    "parent_id": active.index.values,
                    "frozen": False,
                    "depth": active.depth.values,
                    "radius": active.radius.values * self.config.radius_taper,
                    "vessel_type": active.vessel_type.values,
                    "layer_id": active.layer_id.values,
                },
                index=to_freeze.index,
            )
            self.particles.update_particles(continuations)

        frozen_originals = pd.DataFrame(
            {"frozen": True, "freeze_time": self.clock()}, index=active.index
        )
        self.particles.update_particles(frozen_originals)

        if len(available) < len(active) * 3:
            self.add_particles()


class PathExtinction(Component):
    """Component for controlling extinction of active paths based force.

    The per-tip threshold is the ``particle.extinction_threshold`` value
    pipeline (base ``force_threshold``), so other components can modify it —
    PerfusionDemand raises it for tips in hypoxic tissue.
    """

    CONFIGURATION_DEFAULTS = {
        "path_extinction": {
            "force_threshold": 10.0,  # Force magnitude threshold for extinction
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return ["frozen", "freeze_time", "path_id", "vx", "vy", "vz"]

    def setup(self, builder: Builder) -> None:
        self.config = builder.configuration.path_extinction
        self.force_threshold = float(self.config.force_threshold)
        self.clock = builder.time.clock()
        self.particles = builder.components.get_components_by_type(Particle3D)[0]

        # Get force pipelines
        self.force_magnitude = builder.value.get_value("particle.force.magnitude")
        self.force_x = builder.value.get_value("particle.force.x")
        self.force_y = builder.value.get_value("particle.force.y")
        self.force_z = builder.value.get_value("particle.force.z")
        self.threshold = builder.value.register_value_producer(
            "particle.extinction_threshold",
            source=lambda index: pd.Series(self.force_threshold, index=index),
        )

    def on_time_step(self, event: Event) -> None:
        pop = self.population_view.get(event.index, self.required_attributes)
        active = pop[~pop.frozen & (pop.path_id >= 0)]

        if active.empty:
            return

        force_values = self.force_magnitude(active.index)
        to_freeze = active[force_values > self.threshold(active.index)]

        if not to_freeze.empty:
            updates = pd.DataFrame(
                {
                    "frozen": True,
                    "freeze_time": self.clock(),
                    "path_id": -1,  # Mark as end of path
                },
                index=to_freeze.index,
            )
            self.particles.update_particles(updates)


class PathSplitter(Component):
    """Component for splitting particle paths into two branches.

    When the parent vessel has a caliber (radius > 0), daughter calibers are
    assigned by Murray's law and the branch angles follow from the radii via
    the minimum-work principle; otherwise the configured ``split_angle`` is
    used as before.
    """

    CONFIGURATION_DEFAULTS = {
        "path_splitter": {
            "split_interval": 200,
            "split_angle": 30,
            "split_probability": 0.5,
            "max_depth": 4,
            "murray_exponent": 3.0,  # r_parent^k = r_major^k + r_minor^k
            "flow_asymmetry": 0.15,  # minor daughter flow fraction in [0.5 - this, 0.5]
            "min_radius": 0.002,  # caliber floor (capillary scale)
            # Split probability scales as (min_radius / radius) ** this, so wide
            # trunks run long between branch points while narrow twigs branch at
            # the full split_probability; 0 restores caliber-independent cadence
            "caliber_cadence_exponent": 0.0,
            # Comb-like side branching: parents wider than side_branch_radius
            # split at the full split_probability (no cadence damping) and
            # split asymmetrically -- the trunk continues at nearly its own
            # caliber while the side branch takes the Murray minor caliber for
            # this flow fraction, leaving at a near-perpendicular angle, on a
            # random side. 0 disables the mode (dichotomous splitting only)
            "side_branch_flow": 0.0,
            "side_branch_radius": 0.008,
            # Combs form after the arcades establish (retinal development
            # spreads the trunks from the disc first; secondary branching
            # follows). The default predates any simulation start = always on
            "side_branch_start_time": "2000-01-01",
            # Per-split-round emission probability for side-branching trunks
            # (replaces split_probability and the cadence damping for them);
            # sets the comb-tooth spacing along the arcades
            "side_branch_probability": 0.5,
            # Keep at least this many active tips per tree: when a tree's
            # growth front thins below the floor, frozen vessels re-sprout to
            # top it up (angiogenic sprouting wherever the front dies, not
            # only at total extinction). 1 keeps the legacy behavior of
            # re-sprouting only when a tree has no active tips at all
            "min_active_tips": 1,
            # The tip floor applies only to an established tree (one with a
            # real front to lose): below this many frozen particles, few tips
            # is the natural early state and topping up just crowds the disc
            "resprout_established_size": 800,
            # Crowding gate: skip a split when the tip already has at least
            # max_crowding frozen neighbors within crowding_radius. Daughters
            # born into saturated space are pushed over the extinction
            # threshold almost immediately, and the resulting extinction
            # cascade is what makes marginal seeds collapse; gating keeps the
            # branching at the growth front. 0 disables the gate (legacy)
            "crowding_radius": 0.06,
            "max_crowding": 0,
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return [
            "x",
            "y",
            "z",
            "vx",
            "vy",
            "vz",
            "frozen",
            "freeze_time",
            "unfreeze_time",
            "depth",
            "parent_id",
            "path_id",
            "radius",
            "vessel_type",
            "layer_id",
        ]

    def setup(self, builder: Builder) -> None:
        self.config = builder.configuration.path_splitter
        self.particles_to_add = builder.configuration.population.population_size
        self.step_count = 0
        self.next_path_id = builder.configuration.particles.initial_circle.n_vessels + 1
        self.step_size = builder.configuration.time.step_size
        self.randomness = builder.randomness.get_stream("path_splitter")
        self.clock = builder.time.clock()
        self.side_branch_start = pd.Timestamp(self.config.side_branch_start_time)
        self.simulant_creator = builder.population.get_simulant_creator()
        self.particles = builder.components.get_components_by_type(Particle3D)[0]
        self.freezer = builder.components.get_components_by_type(PathFreezer)[0]

    def add_particles(self):
        self.simulant_creator(self.particles_to_add)

    def on_time_step(self, event: Event) -> None:
        self.step_count += 1
        if self.step_count % self.config.split_interval == 0:
            pop = self.population_view.get(event.index, self.required_attributes)
            self.split_paths(pop)

    def split_paths(self, pop: pd.DataFrame) -> None:
        """Split active tips into branches, and re-sprout trees with no tips.

        Active tips split in place (freezing the original). Separately, any
        vessel type whose tree has no active tip left sprouts new branches off
        its own frozen vessels — without this, a tree whose tips all went
        extinct could never grow again while the other tree keeps splitting.
        """
        updates = []
        active = pop[~pop.frozen & (pop.path_id >= 0)]

        if not active.empty:
            to_consider = self.randomness.filter_for_probability(
                active.index, self.split_probabilities(active)
            )
            to_split = self.eligible(pop, to_consider)
            if not to_split.empty:
                updates.extend(self.split_unfrozen(pop, to_split) or [])

        # Per-type re-sprouting from frozen vessels (angiogenic sprouting):
        # a tree whose growth front thins below min_active_tips tops itself
        # up from its own frozen vessels; a tree with no tips at all
        # re-sprouts unconditionally (the legacy bootstrap)
        frozen_on_path = pop[pop.frozen & (pop.path_id >= 0)]
        for vessel_type in np.unique(frozen_on_path.vessel_type):
            n_active = int((active.vessel_type == vessel_type).sum())
            floor = int(self.config.min_active_tips)
            candidates = frozen_on_path[frozen_on_path.vessel_type == vessel_type]
            # The tip floor only applies to an established tree (one with a
            # real front to lose): early on, few tips is the natural state,
            # and topping up then just crowds the disc
            if len(candidates) < int(self.config.resprout_established_size):
                floor = 1
            if n_active >= floor:
                continue
            to_consider = self.randomness.filter_for_probability(
                candidates.index, 0.01, f"active_empty_{vessel_type}"
            )
            to_split = self.eligible(pop, to_consider)
            if n_active > 0:
                to_split = to_split[: floor - n_active]
            if not to_split.empty:
                updates.extend(self.split_frozen(pop, to_split) or [])

        self.commit(updates)

    def eligible(self, pop: pd.DataFrame, candidates: pd.Index) -> pd.Index:
        """Split candidates below the depth ceiling and past the crowding gate."""
        not_too_deep = pop.loc[candidates, "depth"].to_numpy() < self.config.max_depth
        return self.uncrowded(pop, candidates[not_too_deep])

    def commit(self, updates: list) -> None:
        """Write a batch of split updates through Particle3D."""
        if updates:
            self.particles.update_particles(pd.concat(updates, axis=0))

    def side_branch_flow_now(self) -> float:
        """The comb mode's flow fraction, zero before its start time."""
        if self.clock() < self.side_branch_start:
            return 0.0
        return float(self.config.side_branch_flow)

    def uncrowded(self, pop: pd.DataFrame, to_split: pd.Index) -> pd.Index:
        """Drop split candidates whose surroundings are already saturated.

        A daughter born where many frozen vessels stand is pushed over the
        extinction threshold almost immediately, and the cascade of such
        births and deaths is what makes marginal seeds collapse. Gating the
        split keeps branching at the growth front, where there is room.
        """
        limit = int(self.config.max_crowding)
        if limit <= 0 or to_split.empty:
            return to_split
        candidates = pop.loc[to_split]
        # Comb emission by wide trunks is exempt: a trunk's own frozen trail
        # dominates its neighbor count, and suppressing teeth where teeth
        # belong is exactly the wrong response to crowding
        if "radius" in candidates.columns and self.side_branch_flow_now() > 0:
            gated = candidates.radius.to_numpy(float) <= float(self.config.side_branch_radius)
        else:
            gated = np.ones(len(candidates), dtype=bool)
        neighbor_lists = self.freezer.query_radius(
            candidates, float(self.config.crowding_radius)
        )
        if neighbor_lists is None:
            return to_split
        counts = np.array([len(neighbors) for neighbors in neighbor_lists])
        return to_split[~gated | (counts < limit)]

    def resprout_at(self, pop: pd.DataFrame, to_split: pd.Index) -> None:
        """Sprout new branches off the given frozen particles, now.

        The developmental wave calls this when its front stalls: sprouts are
        aimed at specific frozen vessels beside unserved tissue, unlike the
        tip floor in split_paths, which tops a tree up anywhere (and, the
        ninth pass found, seeds cascades in healthy regions doing so).
        Respects the depth ceiling and the crowding gate.
        """
        to_split = self.eligible(pop, to_split)
        if not to_split.empty:
            self.commit(self.split_frozen(pop, to_split) or [])

    def split_probabilities(self, active: pd.DataFrame) -> pd.Series:
        """Per-tip split probability, reduced for wide-caliber tips.

        Real vessel segments keep a roughly constant length-to-diameter ratio,
        so trunks run long between branch points while capillaries branch
        densely. With ``caliber_cadence_exponent`` > 0, the configured
        ``split_probability`` applies at the ``min_radius`` caliber floor and
        wider tips split less often, scaled by
        ``(min_radius / radius) ** caliber_cadence_exponent``. Uncalibered
        tips (radius <= 0) keep the base probability.
        """
        base = float(self.config.split_probability)
        exponent = float(self.config.caliber_cadence_exponent)
        radii = active.radius.to_numpy(dtype=float)
        factors = np.ones(len(radii))
        if exponent != 0.0:
            calibered = radii > 0
            factors[calibered] = (
                float(self.config.min_radius) / radii[calibered]
            ) ** exponent
        probabilities = np.clip(base * factors, 0.0, 1.0)
        # Side-branching trunks are exempt from the cadence damping: real
        # arcades emit side branches at short, comb-like intervals, at
        # their own emission probability
        if self.side_branch_flow_now() > 0:
            wide = radii > float(self.config.side_branch_radius)
            probabilities[wide] = float(self.config.side_branch_probability)
        return pd.Series(probabilities, index=active.index)

    def split_frozen(self, pop, to_split):
        available = pop[~pop.frozen & (pop.path_id < 0)]
        if len(available) < len(to_split):
            self.add_particles()
            return

        new_branches = available.iloc[: len(to_split)]

        angle_rad = np.radians(90) * self.randomness.choice(
            to_split, [-1, 1], [0.5, 0.5], "split_direction"
        )
        angle_rad = angle_rad * (
            0.75 + 0.5 * (self.randomness.get_draw(to_split, "split_angle"))
        )

        # A side branch off a frozen trunk takes the Murray minor-daughter caliber
        _, side_radii, _, _ = self.split_radii_and_angles(pop, to_split)

        # Track updates for frozen originals and new branches
        updates = []

        for idx, orig_idx in enumerate(to_split):
            original = pop.loc[orig_idx]
            vel = np.array([original.vx, original.vy, original.vz])
            speed = np.linalg.norm(vel)
            if speed == 0:
                continue

            # Calculate normalized velocity and perpendicular vector
            vel_norm = vel / speed
            perp = np.array([0, -vel_norm[2], vel_norm[1]])
            if np.allclose(perp, 0):
                perp = np.array([-vel_norm[1], vel_norm[0], 0])
            perp = perp / np.linalg.norm(perp)

            # Calculate new velocities for both branches
            rot_matrix_1 = self._rotation_matrix(perp, angle_rad[orig_idx])
            new_vel_1 = rot_matrix_1 @ vel

            # Normalize new velocities for position offsets
            new_vel_1_norm = new_vel_1 / np.linalg.norm(new_vel_1)

            # Calculate offset positions
            original_pos = np.array([original.x, original.y, original.z])
            pos_1 = original_pos + new_vel_1_norm * speed * self.step_size

            new_branch_1 = pd.DataFrame(
                {
                    "x": [pos_1[0]],
                    "y": [pos_1[1]],
                    "z": [pos_1[2]],
                    "vx": [new_vel_1[0]],
                    "vy": [new_vel_1[1]],
                    "vz": [new_vel_1[2]],
                    "frozen": [False],
                    "freeze_time": [pd.NaT],
                    "depth": [original.depth + 1],
                    "path_id": [self.next_path_id],
                    "parent_id": [orig_idx],
                    "radius": [side_radii[orig_idx]],
                    "vessel_type": [original.vessel_type],
                    "layer_id": [original.layer_id],
                },
                index=[new_branches.iloc[idx].name],
            )
            updates.append(new_branch_1)

        return updates

    def split_radii_and_angles(self, pop, to_split):
        """Murray-law daughter radii and branch angles for each split point.

        Splits with an uncalibered parent (radius <= 0) fall back to the
        configured ``split_angle`` and inherit zero radius.
        """
        parent_radii = pop.loc[to_split, "radius"].astype(float)
        exponent = float(self.config.murray_exponent)
        min_radius = float(self.config.min_radius)

        # Minor daughter carries flow fraction f in [0.5 - flow_asymmetry, 0.5]
        draw = self.randomness.get_draw(to_split, "flow_fraction")
        flow_fractions = 0.5 - self.config.flow_asymmetry * draw
        side_flow = self.side_branch_flow_now()
        side_branching = pd.Series(False, index=to_split)
        if side_flow > 0:
            # Comb mode: a wide trunk keeps nearly its own caliber and sheds
            # a small side branch (flow fraction ~ side_branch_flow, +/-25%)
            side_branching = parent_radii > float(self.config.side_branch_radius)
            flow_fractions = flow_fractions.where(
                ~side_branching, side_flow * (0.75 + 0.5 * draw)
            )
        major_radii, minor_radii = murray_daughter_radii(
            parent_radii, flow_fractions, exponent
        )
        major_radii = major_radii.clip(lower=min_radius)
        minor_radii = minor_radii.clip(lower=min_radius)

        with np.errstate(invalid="ignore", divide="ignore"):
            murray_major, murray_minor = murray_bifurcation_angles(
                parent_radii.clip(lower=min_radius), major_radii, minor_radii
            )

        # Fall back to the configured split angle where the parent has no caliber
        fallback = np.radians(self.config.split_angle / 2)
        has_caliber = parent_radii.values > 0
        angle_major = pd.Series(np.where(has_caliber, murray_major, fallback), index=to_split)
        angle_minor = pd.Series(np.where(has_caliber, murray_minor, fallback), index=to_split)
        major_radii = major_radii.where(has_caliber, 0.0)
        minor_radii = minor_radii.where(has_caliber, 0.0)

        # Retain some stochasticity around the optimal angles
        noise_1 = 0.75 + 0.5 * self.randomness.get_draw(to_split, "split_angle")
        noise_2 = 0.75 + 0.5 * self.randomness.get_draw(to_split, "split_angle_2")
        angle_major = angle_major * noise_1
        angle_minor = angle_minor * noise_2
        if side_flow > 0 and bool(side_branching.any()):
            # Comb teeth leave on a random side; mirroring both angles keeps
            # the trunk's small deviation opposite the tooth
            signs = self.randomness.choice(
                to_split, [-1.0, 1.0], [0.5, 0.5], "side_branch_side"
            )
            signs = signs.astype(float).where(side_branching, 1.0)
            angle_major = angle_major * signs
            angle_minor = angle_minor * signs
        return major_radii, minor_radii, angle_major, angle_minor

    def split_unfrozen(self, pop, to_split):
        available = pop[~pop.frozen & (pop.path_id < 0)]
        if len(available) < 2 * len(to_split):
            self.add_particles()
            return

        # Sample particles for new branches - two per split point
        new_branches = available.iloc[: (2 * len(to_split))]
        major_radii, minor_radii, angle_rad_1, angle_rad_2 = self.split_radii_and_angles(
            pop, to_split
        )

        # Track updates for frozen originals and new branches
        updates = []

        for idx, orig_idx in enumerate(to_split):
            original = pop.loc[orig_idx]
            vel = np.array([original.vx, original.vy, original.vz])
            speed = np.linalg.norm(vel)
            if speed == 0:
                continue

            # Calculate normalized velocity and perpendicular vector
            vel_norm = vel / speed
            perp = np.array([0, -vel_norm[2], vel_norm[1]])
            if np.allclose(perp, 0):
                perp = np.array([-vel_norm[1], vel_norm[0], 0])
            perp = perp / np.linalg.norm(perp)

            # Calculate new velocities for both branches
            rot_matrix_1 = self._rotation_matrix(perp, angle_rad_1[orig_idx])
            rot_matrix_2 = self._rotation_matrix(perp, -angle_rad_2[orig_idx])
            new_vel_1 = rot_matrix_1 @ vel
            new_vel_2 = rot_matrix_2 @ vel

            # Normalize new velocities for position offsets
            new_vel_1_norm = new_vel_1 / np.linalg.norm(new_vel_1)
            new_vel_2_norm = new_vel_2 / np.linalg.norm(new_vel_2)

            # Calculate offset positions
            original_pos = np.array([original.x, original.y, original.z])
            pos_1 = original_pos + new_vel_1_norm * speed * self.step_size
            pos_2 = original_pos + new_vel_2_norm * speed * self.step_size

            # Create DataFrame rows with correct dtypes from the start
            # Freeze original particle at split point
            original_update = pd.DataFrame(
                {
                    "x": [original.x],
                    "y": [original.y],
                    "z": [original.z],
                    "vx": [original.vx],
                    "vy": [original.vy],
                    "vz": [original.vz],
                    "frozen": [True],
                    "freeze_time": [self.clock()],
                    "depth": [original.depth],
                    "path_id": [original.path_id],
                    "parent_id": [original.parent_id],
                    "radius": [original.radius],
                    "vessel_type": [original.vessel_type],
                    "layer_id": [original.layer_id],
                },
                index=[orig_idx],
            )
            updates.append(original_update)

            new_branch_1 = pd.DataFrame(
                {
                    "x": [pos_1[0]],
                    "y": [pos_1[1]],
                    "z": [pos_1[2]],
                    "vx": [new_vel_1[0]],
                    "vy": [new_vel_1[1]],
                    "vz": [new_vel_1[2]],
                    "frozen": [False],
                    "freeze_time": [pd.NaT],
                    "depth": [original.depth],
                    "path_id": [self.next_path_id],
                    "parent_id": [orig_idx],
                    "radius": [major_radii[orig_idx]],
                    "vessel_type": [original.vessel_type],
                    "layer_id": [original.layer_id],
                },
                index=[new_branches.iloc[2 * idx].name],
            )
            updates.append(new_branch_1)

            new_branch_2 = pd.DataFrame(
                {
                    "x": [pos_2[0]],
                    "y": [pos_2[1]],
                    "z": [pos_2[2]],
                    "vx": [new_vel_2[0]],
                    "vy": [new_vel_2[1]],
                    "vz": [new_vel_2[2]],
                    "frozen": [False],
                    "freeze_time": [pd.NaT],
                    "depth": [original.depth + 1],
                    "path_id": [self.next_path_id],
                    "parent_id": [orig_idx],
                    "radius": [minor_radii[orig_idx]],
                    "vessel_type": [original.vessel_type],
                    "layer_id": [original.layer_id],
                },
                index=[new_branches.iloc[2 * idx + 1].name],
            )
            updates.append(new_branch_2)

            self.next_path_id += 2
        return updates

    @staticmethod
    def _rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
        """Return the rotation matrix for rotation around axis by theta radians."""
        axis = axis / np.sqrt(np.dot(axis, axis))
        a = np.cos(theta / 2.0)
        b, c, d = -axis * np.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array(
            [
                [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
            ]
        )


def anastomosis_targets(
    tips: pd.DataFrame,
    frozen: pd.DataFrame,
    neighbor_lists,
    max_target_radius: float,
    min_layer: int = 0,
) -> pd.Series:
    """Match each tip to the nearest opposite-tree capillary within reach.

    ``neighbor_lists[i]`` holds positional indices into ``frozen`` for
    ``tips.iloc[i]`` (as returned by PathFreezer's KDTree query). A frozen
    particle qualifies as a target when it belongs to the *other* tree
    (both vessel types positive and different) and its caliber is at most
    ``max_target_radius`` — capillaries join capillaries, tips don't fuse
    into trunks. Returns a Series mapping tip index to target index for
    the tips that found a target.
    """
    matches = {}
    tip_positions = tips[["x", "y", "z"]].to_numpy(dtype=float)
    for i, neighbors in enumerate(neighbor_lists):
        if len(neighbors) == 0:
            continue
        candidates = frozen.iloc[list(neighbors)]
        tip_type = tips.vessel_type.iloc[i]
        candidates = candidates[
            (candidates.vessel_type > 0)
            & (candidates.vessel_type != tip_type)
            & (candidates.radius > 0)
            & (candidates.radius <= max_target_radius)
        ]
        if min_layer > 0:
            candidates = candidates[candidates.layer_id >= min_layer]
        if candidates.empty:
            continue
        offsets = candidates[["x", "y", "z"]].to_numpy(dtype=float) - tip_positions[i]
        distances = np.linalg.norm(offsets, axis=1)
        matches[tips.index[i]] = candidates.index[int(np.argmin(distances))]
    return pd.Series(matches, dtype=int)


class PathAnastomosis(Component):
    """Fuses capillary-caliber tips onto the other tree, closing the loops.

    Trees don't perfuse; circuits do (roadmap idea 4). When an active
    capillary-caliber tip comes within ``capture_radius`` of the *other*
    tree's capillary-caliber frozen segments, the tip freezes and records
    the join in the ``anastomosis_id`` column, turning the two trees into
    one perfusable graph. Uses PathFreezer's KDTree, so eligible targets
    refresh at the freezer's cadence.
    """

    CONFIGURATION_DEFAULTS = {
        "path_anastomosis": {
            "capture_radius": 0.03,  # tip-to-target distance that triggers fusion
            "max_tip_radius": 0.004,  # only capillary-caliber tips anastomose
            "max_target_radius": 0.004,  # and only onto capillary-caliber segments
            "probability": 0.5,  # per-step fusion probability once in range
            # Fuse only in plexus layers at or below this index. Real capillary
            # loop closure lives in the deeper plexuses; superficial vessels at
            # fundus resolution read as trees. 0 fuses anywhere (legacy)
            "min_layer": 0,
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return [
            "x",
            "y",
            "z",
            "frozen",
            "freeze_time",
            "path_id",
            "radius",
            "vessel_type",
            "anastomosis_id",
            "layer_id",
        ]

    def setup(self, builder: Builder) -> None:
        self.config = builder.configuration.path_anastomosis
        self.clock = builder.time.clock()
        self.randomness = builder.randomness.get_stream("path_anastomosis")
        self.freezer = builder.components.get_components_by_type(PathFreezer)[0]
        self.particles = builder.components.get_components_by_type(Particle3D)[0]

    def on_time_step(self, event: Event) -> None:
        pop = self.population_view.get(event.index, self.required_attributes)
        tips = pop[
            ~pop.frozen
            & (pop.path_id >= 0)
            & (pop.vessel_type > 0)
            & (pop.radius > 0)
            & (pop.radius <= self.config.max_tip_radius)
            & (pop.layer_id >= int(self.config.min_layer))
        ]
        if tips.empty:
            return

        neighbor_lists = self.freezer.query_radius(tips, float(self.config.capture_radius))
        if neighbor_lists is None:
            return

        frozen = self.freezer.frozen_particles()
        targets = anastomosis_targets(
            tips,
            frozen,
            neighbor_lists,
            float(self.config.max_target_radius),
            int(self.config.min_layer),
        )
        if targets.empty:
            return

        to_join = self.randomness.filter_for_probability(
            targets.index, self.config.probability
        )
        if to_join.empty:
            return

        updates = pd.DataFrame(
            {
                "frozen": True,
                "freeze_time": self.clock(),
                "path_id": -1,  # the path ends here, like extinction
                "anastomosis_id": targets[to_join],
            },
            index=to_join,
        )
        self.particles.update_particles(updates)


class PathDLA(Component):
    """Component for freezing particles at the end of a path using DLA.

    The near radius scales exponentially from initial_near_radius to final_near_radius
    between dla_start_time and dla_end_time.
    """

    CONFIGURATION_DEFAULTS = {
        "path_dla": {
            "stickiness": 0.9,
            "initial_near_radius": 0.1,
            "final_near_radius": 0.01,
            "dla_start_time": "2000-01-01",  # Start time for DLA freezing
            "dla_end_time": "2001-01-01",  # End time for radius scaling
            "attach_radius": 0.002,  # caliber assigned to DLA-attached particles
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return [
            "x",
            "y",
            "z",
            "frozen",
            "freeze_time",
            "depth",
            "path_id",
            "parent_id",
            "vessel_type",
            "layer_id",
        ]

    def setup(self, builder: Builder) -> None:
        """Setup the component with configuration and validate parameters."""
        self.config = builder.configuration.path_dla
        self.randomness = builder.randomness.get_stream("path_dla")
        self.clock = builder.time.clock()
        self.freezer = builder.components.get_components_by_type(PathFreezer)[0]
        self.particles = builder.components.get_components_by_type(Particle3D)[0]

        # Convert times to pandas Timestamps
        self.dla_start_time = pd.Timestamp(self.config.dla_start_time)
        self.dla_end_time = pd.Timestamp(self.config.dla_end_time)

        # Validate configuration
        if self.dla_end_time <= self.dla_start_time:
            raise ValueError("dla_end_time must be after dla_start_time")

        if self.config.initial_near_radius <= 0 or self.config.final_near_radius <= 0:
            raise ValueError("near radius values must be positive")

        if self.config.final_near_radius > self.config.initial_near_radius:
            raise ValueError("final_near_radius must be smaller than initial_near_radius")

        # Calculate decay rate for exponential scaling
        total_time = (self.dla_end_time - self.dla_start_time).total_seconds()
        self.decay_rate = (
            -np.log(self.config.final_near_radius / self.config.initial_near_radius)
            / total_time
        )

    def get_current_near_radius(self) -> float:
        """Calculate the current near radius based on exponential decay."""
        current_time = self.clock()

        if current_time < self.dla_start_time:
            return self.config.initial_near_radius
        elif current_time > self.dla_end_time:
            return self.config.final_near_radius

        # Calculate time since start
        time_elapsed = (current_time - self.dla_start_time).total_seconds()

        # Calculate exponentially decayed radius
        current_radius = self.config.initial_near_radius * np.exp(
            -self.decay_rate * time_elapsed
        )
        return current_radius

    def on_time_step(self, event: Event) -> None:
        """Perform DLA freezing with current near radius if after start time."""
        if self.clock() >= self.dla_start_time:
            self.near_radius = self.get_current_near_radius()
            pop = self.population_view.get(event.index, self.required_attributes)
            self.dla_freeze(pop)

    def update_tree(self, pop):
        self._current_frozen = pop[pop.frozen]
        if len(self._current_frozen) < 2:
            self._current_tree = None
        else:
            self._current_tree = cKDTree(self._current_frozen[["x", "y", "z"]].values)

    def dla_freeze(self, pop: pd.DataFrame) -> None:
        """Freeze particles near frozen particles using DLA.
        Only freeze to particles with path_id < 0
        """
        #  only use particles with path_id < 0 (i.e. in frozen DataFrame, not all in freezer object )
        frozen = pop[pop.frozen]  # & (pop.path_id < 0)]
        if frozen.empty:
            return
        self.update_tree(frozen)

        not_frozen = pop[~pop.frozen & (pop.path_id < 0)]
        if not_frozen.empty:
            return

        near_frozen_indices = self._current_tree.query_ball_point(
            not_frozen[["x", "y", "z"]].values, self.near_radius
        )
        near_particles = np.array([len(indices) > 0 for indices in near_frozen_indices])
        stickiness_probabilities = self.randomness.get_draw(
            not_frozen.index, additional_key="stickiness"
        )

        freeze_condition = stickiness_probabilities < self.config.stickiness
        freeze_mask = near_particles & freeze_condition

        to_freeze = not_frozen[freeze_mask]
        if not to_freeze.empty:
            attachment_positions = [
                indices[0] for indices in near_frozen_indices[freeze_mask]
            ]
            updates = pd.DataFrame(
                {
                    "parent_id": frozen.index[attachment_positions],
                    "path_id": 1,
                    "depth": 1000,
                    "frozen": False,
                    "freeze_time": pd.NaT,
                    "radius": self.config.attach_radius,
                    "vessel_type": frozen.vessel_type.values[attachment_positions],
                    "layer_id": frozen.layer_id.values[attachment_positions],
                },
                index=to_freeze.index,
            )
            self.particles.update_particles(updates)

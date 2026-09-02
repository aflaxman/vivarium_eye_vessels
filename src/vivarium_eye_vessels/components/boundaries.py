from typing import Dict, List, Protocol

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from vivarium import Component
from vivarium.framework.engine import Builder
from vivarium.framework.event import Event

from vivarium_eye_vessels.components.particles import (
    VESSEL_TYPE_ARTERY,
    VESSEL_TYPE_VEIN,
    PathFreezer,
    PathSplitter,
)


class ForceCalculator(Protocol):
    """Protocol defining the interface for force calculation strategies"""

    def calculate_force_magnitude(self, distances: np.ndarray) -> np.ndarray:
        """Calculate force magnitudes based on distances"""
        pass


class HookeanForce:
    """Implements Hooke's law force calculation"""

    def __init__(self, spring_constant: float):
        self.spring_constant = spring_constant

    def calculate_force_magnitude(self, distances: np.ndarray) -> np.ndarray:
        return self.spring_constant * distances


class MagneticForce:
    """Implements inverse square law force calculation"""

    def __init__(self, magnetic_strength: float, min_distance: float):
        self.magnetic_strength = magnetic_strength
        self.min_distance = min_distance

    def calculate_force_magnitude(self, distances: np.ndarray) -> np.ndarray:
        capped_distances = np.maximum(distances, self.min_distance)
        return self.magnetic_strength / (capped_distances * capped_distances)


class BaseForceComponent(Component):
    """Base class for force-based components with shared caching logic"""

    @property
    def required_attributes(self) -> List[str]:
        return ["x", "y", "z", "frozen"]

    @property
    def filter_str(self) -> str:
        return "not frozen"

    def setup(self, builder: Builder) -> None:
        self.force_cache = {}
        self.clock = builder.time.clock()

        # Register force modifiers
        for axis in ["x", "y", "z"]:
            builder.value.register_value_modifier(
                f"particle.force.{axis}",
                modifier=getattr(self, f"force_{axis}"),
                required_resources=self.required_attributes,
            )

    def setup_force_calculator(self, config: Dict) -> None:
        # Set up force calculator
        if config.force_type == "magnetic":
            self.force_calculator = MagneticForce(
                float(config.magnetic_strength), float(config.min_distance)
            )
        else:
            self.force_calculator = HookeanForce(float(config.spring_constant))

    def get_cached_forces(self, index: pd.Index) -> np.ndarray:
        """Get cached forces or calculate them if needed"""
        current_time = self.clock()
        cache_key = (current_time, tuple(index))

        if cache_key not in self.force_cache:
            pop = self.population_view.get(index, self.required_attributes)
            active_particles = pop.query(self.filter_str)

            if active_particles.empty:
                self.force_cache[cache_key] = np.zeros((len(index), 3))
            else:
                forces = np.zeros((len(index), 3))
                active_forces = self.calculate_forces_vectorized(active_particles)
                forces[index.get_indexer(active_particles.index)] = active_forces
                self.force_cache[cache_key] = forces

            # Clear old cache entries
            self.force_cache = {
                k: v for k, v in self.force_cache.items() if k[0] == current_time
            }

        return self.force_cache[cache_key]

    def force_x(self, index: pd.Index, forces: pd.Series) -> pd.Series:
        forces += pd.Series(self.get_cached_forces(index)[:, 0], index=index)
        return forces

    def force_y(self, index: pd.Index, forces: pd.Series) -> pd.Series:
        forces += pd.Series(self.get_cached_forces(index)[:, 1], index=index)
        return forces

    def force_z(self, index: pd.Index, forces: pd.Series) -> pd.Series:
        forces += pd.Series(self.get_cached_forces(index)[:, 2], index=index)
        return forces


class EllipsoidContainment(BaseForceComponent):
    """Component that keeps particles within an ellipsoid boundary"""

    CONFIGURATION_DEFAULTS = {
        "ellipsoid_containment": {
            "a": 1.0,
            "b": 1.0,
            "c": 1.0,
            "force_type": "hookean",  # "magnetic" or "hookean"
            "magnetic_strength": 0.1,
            "min_distance": 0.01,
            "spring_constant": 1.0,
        }
    }

    def setup(self, builder: Builder) -> None:
        super().setup(builder)
        config = builder.configuration.ellipsoid_containment
        self.setup_force_calculator(config)

        # Set up geometry parameters
        self.a = float(config.a)
        self.b = float(config.b)
        self.c = float(config.c)
        self.a2 = self.a * self.a
        self.b2 = self.b * self.b
        self.c2 = self.c * self.c

    def calculate_forces_vectorized(self, particles: pd.DataFrame) -> np.ndarray:
        positions = particles[["x", "y", "z"]].to_numpy()

        # Calculate normalized coordinates
        x_norm = positions[:, 0] / self.a
        y_norm = positions[:, 1] / self.b
        z_norm = positions[:, 2] / self.c

        # Calculate ellipsoid equation value
        ellipsoid_val = x_norm**2 + y_norm**2 + z_norm**2

        # Calculate gradient components
        grad = np.column_stack(
            [2 * x_norm / self.a, 2 * y_norm / self.b, 2 * z_norm / self.c]
        )

        # Initialize forces array
        forces = np.zeros_like(positions)
        outside_mask = ellipsoid_val > 1

        if np.any(outside_mask):
            grad_outside = grad[outside_mask]
            grad_norms = np.linalg.norm(grad_outside, axis=1, keepdims=True)
            normalized_grads = grad_outside / grad_norms

            # Calculate distances from surface
            distances = np.sqrt(ellipsoid_val[outside_mask]) - 1

            # Calculate force magnitudes using the selected force calculator
            force_magnitudes = self.force_calculator.calculate_force_magnitude(distances)

            # Calculate final forces
            forces[outside_mask] = -normalized_grads * force_magnitudes[:, np.newaxis]

        return forces


class CylinderExclusion(BaseForceComponent):
    """Component that repels particles from inside a cylindrical exclusion zone"""

    CONFIGURATION_DEFAULTS = {
        "cylinder_exclusion": {
            "radius": 1.0,
            "center": [0.0, 0.0, 0.0],
            "direction": [0.0, 0.0, 1.0],
            "force_type": "hookean",  # or "magnetic"
            "magnetic_strength": 0.1,
            "min_distance": 0.01,
            "spring_constant": 0.1,
        }
    }

    def setup(self, builder: Builder) -> None:
        super().setup(builder)

        config = builder.configuration.cylinder_exclusion
        self.setup_force_calculator(config)

        # Set up geometry parameters
        self.radius = float(config.radius)
        self.center = np.array(config.center, dtype=float)
        self.direction = np.array(config.direction, dtype=float)
        self.direction /= np.linalg.norm(self.direction)

        # Pre-compute random perpendicular vector
        random_perpendicular = (
            np.array([1, 0, 0]) if abs(self.direction[0]) < 0.9 else np.array([0, 1, 0])
        )
        self.default_outward = np.cross(self.direction, random_perpendicular)
        self.default_outward /= np.linalg.norm(self.default_outward)

    def calculate_forces_vectorized(self, particles: pd.DataFrame) -> np.ndarray:
        positions = particles[["x", "y", "z"]].to_numpy()

        # Calculate relative positions and components
        rel_positions = positions - self.center
        axial_dots = np.dot(rel_positions, self.direction)
        axial_components = axial_dots[:, np.newaxis] * self.direction
        radial_vectors = rel_positions - axial_components
        radial_distances = np.linalg.norm(radial_vectors, axis=1)

        # Calculate penetration depths
        penetrations = self.radius - radial_distances

        # Handle points on axis
        mask_on_axis = radial_distances < 1e-10
        outward_directions = np.zeros_like(positions)
        outward_directions[mask_on_axis] = self.default_outward

        # Calculate outward directions for off-axis points
        mask_off_axis = ~mask_on_axis
        outward_directions[mask_off_axis] = (
            radial_vectors[mask_off_axis] / radial_distances[mask_off_axis, np.newaxis]
        )

        # Apply forces only inside cylinder
        mask_inside = penetrations > 0
        force_magnitudes = np.zeros_like(radial_distances)
        force_magnitudes[mask_inside] = self.force_calculator.calculate_force_magnitude(
            penetrations[mask_inside]
        )

        return outward_directions * force_magnitudes[:, np.newaxis]


class PointRepulsion(BaseForceComponent):
    """Component that creates a point-based repulsion force"""

    CONFIGURATION_DEFAULTS = {
        "point_repulsion": {
            "position": {"x": 0.0, "y": 0.0, "z": 0.0},
            "force_type": "magnetic",  # or "hookean"
            "magnetic_strength": 0.05,
            "min_distance": 0.1,
            "spring_constant": 0.1,
            "radius": 0.05,  # Interaction radius
        }
    }

    def setup(self, builder: Builder) -> None:
        super().setup(builder)

        config = builder.configuration.point_repulsion
        self.setup_force_calculator(config)

        self.position = np.array(
            [float(config.position.x), float(config.position.y), float(config.position.z)]
        )
        self.radius = float(config.radius)

    def calculate_forces_vectorized(self, particles: pd.DataFrame) -> np.ndarray:
        positions = particles[["x", "y", "z"]].to_numpy()

        # Calculate displacements and distances
        displacements = self.position - positions
        distances = np.sqrt(np.sum(displacements**2, axis=1))
        distances = np.where(distances > self.radius, 0, distances)

        # Calculate normalized directions
        with np.errstate(invalid="ignore", divide="ignore"):
            direction_vectors = displacements / distances[:, np.newaxis]
        direction_vectors = np.nan_to_num(direction_vectors)

        # Calculate force magnitudes using selected force calculator
        force_magnitudes = self.force_calculator.calculate_force_magnitude(distances)

        # Return repulsive forces
        return -direction_vectors * force_magnitudes[:, np.newaxis]


class FrozenRepulsion(BaseForceComponent):
    """Component that repels active particles from frozen particles using spatial indexing.

    When particles carry artery/vein identities, repulsion from the *other*
    tree's frozen vessels is scaled by ``cross_type_factor``. A factor below 1
    lets arteries and veins tolerate each other's proximity while still
    avoiding their own tree, which produces the interdigitating arcade
    pattern of the retina.
    """

    CONFIGURATION_DEFAULTS = {
        "frozen_repulsion": {
            "interaction_radius": 0.2,
            "freeze_radius": 0.05,
            "force_type": "magnetic",  # or "hookean"
            "magnetic_strength": 0.1,
            "min_distance": 0.01,
            "spring_constant": 0.1,
            "delay": 1.0,  # days frozen before exerting force on particles in same path
            "cross_type_factor": 1.0,  # repulsion multiplier between artery and vein
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return super().required_attributes + [
            "freeze_time",
            "path_id",
            "parent_id",
            "vessel_type",
        ]

    @property
    def filter_str(self) -> str:
        return "not frozen and path_id >= 0"

    def setup(self, builder: Builder) -> None:
        super().setup(builder)
        config = builder.configuration.frozen_repulsion
        self.setup_force_calculator(config)
        self.clock = builder.time.clock()

        self.interaction_radius = float(config.interaction_radius)
        self.freeze_radius = float(config.freeze_radius)
        self.delay = float(config.delay)
        self.cross_type_factor = float(config.cross_type_factor)
        self.freezer = builder.components.get_components_by_type(PathFreezer)[0]

    def calculate_forces_vectorized(self, particles: pd.DataFrame) -> np.ndarray:
        """Calculate repulsion forces from frozen particles"""
        positions = particles[["x", "y", "z"]].to_numpy()

        forces = np.zeros_like(positions)
        neighbor_lists = self.freezer.query_radius(positions, self.interaction_radius)

        if neighbor_lists is None:
            return forces

        for i, frozen_neighbors in enumerate(neighbor_lists):
            # Calculate displacement vectors from frozen particles
            frozen = self.freezer.get_population(frozen_neighbors)
            frozen = frozen[
                (
                    (frozen.path_id == particles.iloc[i].path_id)
                    & (
                        (self.clock() - frozen.freeze_time) / pd.Timedelta(days=1)
                        > self.delay
                    )
                )
                | (frozen.path_id != particles.iloc[i].path_id)
            ]

            frozen_neighbor_positions = frozen[["x", "y", "z"]].to_numpy()
            displacements = positions[i] - frozen_neighbor_positions

            # Calculate distances
            distances = np.sqrt(np.sum(displacements**2, axis=1))

            # Calculate normalized direction vectors
            with np.errstate(invalid="ignore", divide="ignore"):
                directions = displacements / distances[:, np.newaxis]
            directions = np.nan_to_num(directions)

            # Calculate and sum forces from all frozen neighbors, with weaker
            # repulsion from the other tree (artery vs. vein)
            force_magnitudes = self.force_calculator.calculate_force_magnitude(
                self.interaction_radius - distances
            )
            same_type = frozen.vessel_type.to_numpy() == particles.iloc[i].vessel_type
            type_factors = np.where(same_type, 1.0, self.cross_type_factor)
            forces[i] = np.sum(
                directions * (force_magnitudes * type_factors)[:, np.newaxis], axis=0
            )

        return forces


#######################################
# Hypoxia-driven growth (roadmap 2)   #
#######################################


def generate_demand_sites(semi_axes: np.ndarray, spacing: float) -> np.ndarray:
    """Lattice of tissue demand sites inside the containment ellipsoid.

    Each site represents a patch of tissue that needs perfusion; sites outside
    the ellipsoid are discarded.
    """
    a, b, c = (float(v) for v in semi_axes)
    xs = np.arange(-a, a + spacing / 2, spacing)
    ys = np.arange(-b, b + spacing / 2, spacing)
    zs = np.arange(-c, c + spacing / 2, spacing)
    grid = np.stack(np.meshgrid(xs, ys, zs, indexing="ij"), axis=-1).reshape(-1, 3)
    inside = ((grid[:, 0] / a) ** 2 + (grid[:, 1] / b) ** 2 + (grid[:, 2] / c) ** 2) <= 1.0
    return grid[inside]


def colonization_forces(
    tip_positions: np.ndarray,
    hypoxic_sites: np.ndarray,
    influence_radius: float,
    magnitude: float,
) -> np.ndarray:
    """Space-colonization attraction of growth tips toward hypoxic tissue.

    Following Runions et al. (2005), each hypoxic site recruits only its
    nearest growth tip (within ``influence_radius``); each tip is pulled in
    the normalized mean direction of the sites it won, scaled by
    ``magnitude``. Because a site pulls exactly one tip, tips naturally
    spread apart and fill unperfused territory instead of clustering.
    """
    forces = np.zeros_like(tip_positions, dtype=float)
    if len(tip_positions) == 0 or len(hypoxic_sites) == 0:
        return forces

    tip_tree = cKDTree(tip_positions)
    distances, nearest_tip = tip_tree.query(hypoxic_sites, k=1)
    in_range = (distances <= influence_radius) & (distances > 1e-12)
    if not np.any(in_range):
        return forces

    offsets = hypoxic_sites[in_range] - tip_positions[nearest_tip[in_range]]
    unit_vectors = offsets / distances[in_range, np.newaxis]
    np.add.at(forces, nearest_tip[in_range], unit_vectors)

    norms = np.linalg.norm(forces, axis=1)
    pulled = norms > 1e-12
    forces[pulled] = magnitude * forces[pulled] / norms[pulled, np.newaxis]
    return forces


class PerfusionDemand(BaseForceComponent):
    """Attracts vessel growth tips toward under-perfused tissue.

    A stand-in for hypoxia-driven VEGF signaling: tissue sites farther than
    ``perfusion_radius`` from any frozen vessel are hypoxic, and each hypoxic
    site recruits its nearest active growth tip. Uses PathFreezer's KDTree of
    frozen particles, so perfusion state refreshes at the freezer's cadence.

    When vessels carry artery/vein identities, demand is computed per tree:
    tissue needs both arterial supply and venous drainage, so a site keeps
    recruiting artery tips until an artery is nearby regardless of how well
    it is drained, and vice versa. This keeps the two trees in balance —
    neither can win territory for both.
    """

    CONFIGURATION_DEFAULTS = {
        "perfusion_demand": {
            "site_spacing": 0.1,
            "perfusion_radius": 0.15,
            "influence_radius": 2.0,
            "magnitude": 0.3,
            # Hypoxia chemotaxis acts on capillary sprout tips, not trunks:
            # tips wider than caliber_reference have the attraction attenuated
            # by (reference / radius) ** caliber_exponent. Exponent 0 keeps
            # the legacy caliber-blind attraction
            "caliber_reference": 0.004,
            "caliber_exponent": 0.0,
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return ["x", "y", "z", "frozen", "path_id", "vessel_type", "radius"]

    @property
    def filter_str(self) -> str:
        return "not frozen and path_id >= 0"

    def setup(self, builder: Builder) -> None:
        super().setup(builder)
        config = builder.configuration.perfusion_demand
        self.perfusion_radius = float(config.perfusion_radius)
        self.influence_radius = float(config.influence_radius)
        self.magnitude = float(config.magnitude)
        self.caliber_reference = float(config.caliber_reference)
        self.caliber_exponent = float(config.caliber_exponent)
        self.freezer = builder.components.get_components_by_type(PathFreezer)[0]
        waves = builder.components.get_components_by_type(DevelopmentalWave)
        self.wave = waves[0] if waves else None

        if "ellipsoid_containment" in builder.components.list_components():
            ellipsoid = builder.configuration.ellipsoid_containment
            semi_axes = np.array([ellipsoid.a, ellipsoid.b, ellipsoid.c], dtype=float)
        else:
            semi_axes = np.ones(3)
        self.sites = generate_demand_sites(semi_axes, float(config.site_spacing))

    def hypoxic_sites(
        self, vessel_type: int | None = None, visible_only: bool = True
    ) -> np.ndarray:
        """Demand sites farther than perfusion_radius from any (typed) vessel.

        With ``vessel_type`` given, only frozen vessels of that type count as
        perfusing; with None, any frozen vessel does. When a DevelopmentalWave
        is present, only the sites its front currently exposes are returned,
        unless ``visible_only`` is False (the whole hypoxic field).
        """
        frozen = self.freezer.frozen_particles()
        if frozen is not None and vessel_type is not None:
            frozen = frozen[frozen.vessel_type == vessel_type]
        if frozen is None or frozen.empty:
            hypoxic = self.sites
        else:
            tree = cKDTree(frozen[["x", "y", "z"]].to_numpy(dtype=float))
            distances, _ = tree.query(self.sites, k=1)
            hypoxic = self.sites[distances > self.perfusion_radius]
        if visible_only and self.wave is not None:
            hypoxic = self.wave.visible(hypoxic)
        return hypoxic

    def calculate_forces_vectorized(self, particles: pd.DataFrame) -> np.ndarray:
        tips = particles[["x", "y", "z"]].to_numpy(dtype=float)
        tip_types = particles["vessel_type"].to_numpy()
        radii = particles["radius"].to_numpy(dtype=float)
        # Tips wider than the wave's gate_caliber_max see the whole hypoxic
        # field (arcade trunks follow the astrocyte template out to the
        # periphery, not local VEGF); narrower tips see only what the front
        # exposes. 0 gates every tip
        exempt_over = 0.0
        if self.wave is not None and self.wave.enabled:
            exempt_over = float(self.wave.config.gate_caliber_max)
        forces = np.zeros_like(tips)
        for vessel_type in np.unique(tip_types):
            selected = tip_types == vessel_type
            key = int(vessel_type) if vessel_type > 0 else None
            hypoxic = self.hypoxic_sites(key)
            if exempt_over > 0:
                wide = selected & (radii > exempt_over)
                selected = selected & (radii <= exempt_over)
                if np.any(wide):
                    forces[wide] = colonization_forces(
                        tips[wide],
                        self.hypoxic_sites(key, visible_only=False),
                        self.influence_radius,
                        self.magnitude,
                    )
            forces[selected] = colonization_forces(
                tips[selected], hypoxic, self.influence_radius, self.magnitude
            )
        # Chemotaxis is a capillary-sprout behavior: attenuate the pull on
        # wide tips so trunks and side branches hold their heading
        if self.caliber_exponent > 0:
            radii = particles["radius"].to_numpy(dtype=float)
            attenuation = np.where(
                radii > self.caliber_reference,
                (self.caliber_reference / np.maximum(radii, 1e-12)) ** self.caliber_exponent,
                1.0,
            )
            forces = forces * attenuation[:, np.newaxis]
        return forces


class DevelopmentalWave(Component):
    """Closed-loop growth front: vascularization advances as a wave.

    Real retinal vascularization spreads from the optic disc as a coordinated
    front — hypoxic signaling is strongest just beyond the vascularized zone,
    extinguished behind it, and not yet active far ahead of it. This
    component owns that front as a radius around the disc:

    - PerfusionDemand exposes only the demand sites within
      ``radius + lookahead`` of the disc, so every tip chases the same
      expanding ring instead of an isotropic field;
    - the front advances by ``wave_speed`` per step only while the tissue
      behind it is actually served (perfused fraction inside the front at or
      above ``advance_threshold``, for each vessel type), so a stalling seed
      concentrates demand at the stall instead of failing silently;
    - a front held for ``hold_resprout_steps`` re-sprouts the stalled tree
      from the frozen vessels nearest its unserved sites — growth pressure
      applied exactly where growth is missing, unlike the global tip floor,
      which the ninth pass found seeds extinction cascades in healthy
      regions when it tops a tree up anywhere.

    The control variable is the actual goal (perfusion behind the front)
    rather than a local proxy, which is what distinguishes this from the
    per-tip rules the ninth pass rejected. Disabled, the component is a
    strict no-op and dynamics are unchanged.
    """

    CONFIGURATION_DEFAULTS = {
        "developmental_wave": {
            "enabled": False,
            "start_radius": 0.15,
            "wave_speed": 0.006,  # advance per step while the tissue behind is served
            "lookahead": 0.25,  # demand visible this far beyond the front
            "advance_threshold": 0.85,  # served fraction inside the front to advance
            "hold_resprout_steps": 15,  # held steps before targeted re-sprouting
            "resprout_count": 2,  # sprouts per stalled tree per trigger
            # "combined" advances on any-vessel service of the tissue behind
            # the front; "per_type" requires every tree to serve it before
            # advancing. Combined is the validated default: the artery tree
            # alone plateaus near 0.85 coverage mid-field, so per_type brakes
            # healthy seeds on the artery ceiling while pruning grinds them
            "advance_rule": "combined",
            # Tips wider than this see the whole hypoxic field rather than
            # only what the front exposes (arcade trunks follow the astrocyte
            # template to the periphery, not local VEGF); 0 gates every tip
            "gate_caliber_max": 0.0,
        }
    }

    def setup(self, builder: Builder) -> None:
        self.config = builder.configuration.developmental_wave
        self.enabled = bool(self.config.enabled)
        self.radius = float(self.config.start_radius)
        self.center = np.array(
            builder.configuration.particles.initial_circle.center, dtype=float
        )
        self.hold_steps: Dict[int, int] = {}
        self.freezer = builder.components.get_components_by_type(PathFreezer)[0]
        demands = builder.components.get_components_by_type(PerfusionDemand)
        self.demand = demands[0] if demands else None
        splitters = builder.components.get_components_by_type(PathSplitter)
        self.splitter = splitters[0] if splitters else None
        if self.enabled and self.demand is None:
            raise ValueError("DevelopmentalWave requires a PerfusionDemand component")

    def visible(self, sites: np.ndarray) -> np.ndarray:
        """The demand sites the current front exposes to growth tips."""
        if not self.enabled or len(sites) == 0:
            return sites
        distances = np.linalg.norm(sites - self.center, axis=1)
        return sites[distances <= self.radius + float(self.config.lookahead)]

    def served_fraction(self, vessel_type: int | None) -> float:
        """Fraction of demand sites inside the front perfused by this tree.

        With ``vessel_type`` None, any frozen vessel counts (the combined
        advance rule).
        """
        sites = self.demand.sites
        inside = sites[np.linalg.norm(sites - self.center, axis=1) <= self.radius]
        if len(inside) == 0:
            return 1.0
        frozen = self.freezer.frozen_particles()
        if frozen is not None and vessel_type is not None:
            frozen = frozen[frozen.vessel_type == vessel_type]
        if frozen is None or frozen.empty:
            return 0.0
        tree = cKDTree(frozen[["x", "y", "z"]].to_numpy(dtype=float))
        distances, _ = tree.query(inside, k=1)
        return float((distances <= self.demand.perfusion_radius).mean())

    def on_time_step(self, event: Event) -> None:
        if not self.enabled:
            return
        threshold = float(self.config.advance_threshold)
        if str(self.config.advance_rule) == "combined":
            checks: list[int | None] = [None]
        else:
            checks = [VESSEL_TYPE_ARTERY, VESSEL_TYPE_VEIN]
        stalled = []
        served = []
        for check in checks:
            fraction = self.served_fraction(check)
            served.append(fraction)
            if fraction < threshold:
                held = self.hold_steps.get(check, 0) + 1
                if held >= int(self.config.hold_resprout_steps):
                    stalled.append(check)
                    held = 0
                self.hold_steps[check] = held
            else:
                self.hold_steps[check] = 0
        if min(served) >= threshold:
            self.radius += float(self.config.wave_speed)
        for check in stalled:
            types = (VESSEL_TYPE_ARTERY, VESSEL_TYPE_VEIN) if check is None else (check,)
            for vessel_type in types:
                self.resprout_toward_stall(vessel_type, event)

    def resprout_toward_stall(self, vessel_type: int, event: Event) -> None:
        """Sprout the stalled tree from frozen vessels nearest unserved tissue.

        A hypoxic site is by definition farther than perfusion_radius from
        every frozen vessel of its type, so the nearest frozen candidates are
        the frontier of the stalled region — the wave applies its growth
        pressure there and nowhere else.
        """
        if self.splitter is None:
            return
        stall_sites = self.demand.hypoxic_sites(vessel_type)
        if len(stall_sites) == 0:
            return
        pop = self.population_view.get(event.index, self.splitter.required_attributes)
        candidates = pop[pop.frozen & (pop.path_id >= 0) & (pop.vessel_type == vessel_type)]
        if candidates.empty:
            return
        tree = cKDTree(stall_sites)
        distances, _ = tree.query(candidates[["x", "y", "z"]].to_numpy(dtype=float), k=1)
        order = np.argsort(distances)
        chosen = candidates.index[order[: int(self.config.resprout_count)]]
        self.splitter.resprout_at(pop, pd.Index(chosen))

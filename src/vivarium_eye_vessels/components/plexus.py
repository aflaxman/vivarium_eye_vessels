"""Stratified vascular plexuses connected by diving vessels (roadmap idea 6).

The retina's vasculature is not one sheet: arteries, veins, and the largest
arterioles live in the superficial vascular plexus, while the intermediate
and deep capillary plexuses are capillary-only networks fed by short
"diving" vessels that plunge vertically between layers. This component
gives each vessel particle a home plexus (the ``layer_id`` column), holds
growth tips near their layer's z-plane with a Hookean spring, and lets
capillary-caliber tips occasionally dive one layer deeper — the frozen
trail left during the transit is the diving vessel.
"""

from typing import List

import numpy as np
import pandas as pd
from vivarium import Component
from vivarium.framework.engine import Builder
from vivarium.framework.event import Event

from vivarium_eye_vessels.components.particles import Particle3D


class PlexusLayers(Component):
    """Holds each vessel in its plexus and lets capillary tips dive deeper.

    ``layer_z`` lists the z-plane of each plexus, superficial first; a
    particle's ``layer_id`` indexes into it. A z-only Hookean force pulls
    every active on-path tip toward its own plane, so wide vessels (which
    never dive) stay superficial while diving tips transit downward. Each
    step, an active tip with caliber at most ``dive_radius`` and a layer
    above the deepest dives one layer with probability
    ``dive_probability`` — matching the anatomy, where only capillaries
    populate the intermediate and deep plexuses.
    """

    CONFIGURATION_DEFAULTS = {
        "plexus_layers": {
            "layer_z": [0.04, 0.0, -0.04],  # superficial, intermediate, deep
            "spring_constant": 6.0,
            # Damping on vertical velocity: without it the spring is an
            # undamped oscillator — tips ping-pong through their plane, the
            # z-velocity saturates the terminal-velocity clamp (z counts 10x
            # in the scaled speed norm), and stalled tips go extinct
            "damping": 5.0,
            # Cap on the layer pull; keep below path_extinction.force_threshold
            # (like perfusion_demand.magnitude) so stratification and diving
            # don't push tips into extinction
            "max_force": 0.3,
            "dive_radius": 0.004,  # only capillary-caliber tips dive
            # Tips narrower than this never dive: capillary sprouts belong
            # to their own layer's bed (CapillaryBed). 0 lets every tip dive
            "min_dive_radius": 0.0,
            "dive_probability": 0.02,  # per-step chance an eligible tip dives
            # A diving vessel becomes a capillary when it reaches its plane:
            # tips of the arteriole class in a deeper layer, within
            # plane_tolerance of the plane, take this caliber and from then
            # on obey the CapillaryBed's rules (no splits, fine hypoxia,
            # anastomosis, regression). The intermediate and deep plexuses
            # are capillary-only in the retina; the vertical connector is the
            # trail the tip left on the way down. 0 = legacy (the diver keeps
            # growing as an arteriole-class vessel in the deep layer)
            "arrival_caliber": 0.0,
            "plane_tolerance": 0.01,
            # Judge the dive caliber on each tree's own scale (artery tips are
            # artery_caliber_ratio narrower than vein tips of the same rank):
            # with one absolute dive_radius, more artery than vein tips
            # qualified to dive, and once divers become capillaries the artery
            # tree bled tips into the deep beds and lost 40% of its coverage
            "type_scaled_dive": False,
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return ["z", "vz", "frozen", "path_id", "radius", "layer_id", "vessel_type"]

    def setup(self, builder: Builder) -> None:
        config = builder.configuration.plexus_layers
        self.layer_z = np.asarray(list(config.layer_z), dtype=float)
        self.spring_constant = float(config.spring_constant)
        self.damping = float(config.damping)
        self.max_force = float(config.max_force)
        self.dive_radius = float(config.dive_radius)
        self.min_dive_radius = float(config.min_dive_radius)
        self.dive_probability = float(config.dive_probability)
        self.arrival_caliber = float(config.arrival_caliber)
        self.plane_tolerance = float(config.plane_tolerance)
        self.type_scaled_dive = bool(config.type_scaled_dive)
        self.artery_caliber_ratio = float(
            builder.configuration.particles.artery_caliber_ratio
        )
        self.randomness = builder.randomness.get_stream("plexus_layers")
        self.particles = builder.components.get_components_by_type(Particle3D)[0]

        builder.value.register_value_modifier(
            "particle.force.z",
            modifier=self.layer_force_z,
            required_resources=self.required_attributes,
        )

    def active_tips(self, pop: pd.DataFrame) -> pd.DataFrame:
        return pop[~pop.frozen & (pop.path_id >= 0) & (pop.layer_id >= 0)]

    def layer_force_z(self, index: pd.Index, forces: pd.Series) -> pd.Series:
        """Damped Hookean pull of each active tip toward its own plexus plane."""
        pop = self.population_view.get(index, self.required_attributes)
        tips = self.active_tips(pop)
        if tips.empty:
            return forces
        layers = np.clip(tips.layer_id.to_numpy(int), 0, len(self.layer_z) - 1)
        pull = self.spring_constant * (
            self.layer_z[layers] - tips.z.to_numpy(float)
        ) - self.damping * tips.vz.to_numpy(float)
        forces[tips.index] += np.clip(pull, -self.max_force, self.max_force)
        return forces

    def on_time_step(self, event: Event) -> None:
        """Send an occasional capillary-caliber tip one plexus deeper."""
        pop = self.population_view.get(event.index, self.required_attributes)
        tips = self.active_tips(pop)
        if self.arrival_caliber > 0:
            self.convert_arrivals(tips)
        dive_radius = np.full(len(tips), self.dive_radius)
        if self.type_scaled_dive and "vessel_type" in tips.columns:
            dive_radius[tips.vessel_type.to_numpy() == 1] *= self.artery_caliber_ratio
        eligible = tips[
            (tips.radius > 0)
            & (tips.radius <= dive_radius)
            & (tips.radius >= self.min_dive_radius)
            & (tips.layer_id < len(self.layer_z) - 1)
        ]
        if eligible.empty:
            return
        divers = self.randomness.filter_for_probability(eligible.index, self.dive_probability)
        if divers.empty:
            return
        self.particles.update_particles(
            pd.DataFrame({"layer_id": pop.loc[divers, "layer_id"] + 1}, index=divers)
        )

    def convert_arrivals(self, tips: pd.DataFrame) -> None:
        """Arteriole-class tips that have reached a deeper plane become capillaries."""
        layers = np.clip(tips.layer_id.to_numpy(int), 0, len(self.layer_z) - 1)
        arrived = (
            (tips.layer_id.to_numpy() > 0)
            & (
                tips.radius.to_numpy(float)
                >= max(self.min_dive_radius, self.arrival_caliber * 1.001)
            )
            & (np.abs(tips.z.to_numpy(float) - self.layer_z[layers]) <= self.plane_tolerance)
        )
        if arrived.any():
            self.particles.update_particles(
                pd.DataFrame({"radius": self.arrival_caliber}, index=tips.index[arrived])
            )

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
            "dive_probability": 0.02,  # per-step chance an eligible tip dives
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return ["z", "vz", "frozen", "path_id", "radius", "layer_id"]

    def setup(self, builder: Builder) -> None:
        config = builder.configuration.plexus_layers
        self.layer_z = np.asarray(list(config.layer_z), dtype=float)
        self.spring_constant = float(config.spring_constant)
        self.damping = float(config.damping)
        self.max_force = float(config.max_force)
        self.dive_radius = float(config.dive_radius)
        self.dive_probability = float(config.dive_probability)
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
        eligible = tips[
            (tips.radius > 0)
            & (tips.radius <= self.dive_radius)
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

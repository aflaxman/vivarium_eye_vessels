from typing import List

from vivarium import Component
from vivarium.framework.engine import Builder
from vivarium.framework.event import Event


class SaveParticles(Component):
    @property
    def required_attributes(self) -> List[str]:
        return ["x", "y", "z", "parent_id", "frozen", "depth"]

    def setup(self, builder: Builder) -> None:
        self.seed = builder.configuration.randomness.random_seed

    def on_simulation_end(self, event: Event) -> None:
        pop = self.population_view.get(event.index, self.required_attributes)
        fname = f"{self.seed}.csv.bz2"
        pop.to_csv(fname)

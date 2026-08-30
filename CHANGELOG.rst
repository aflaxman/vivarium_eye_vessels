**v0.1.1 - 08/29/26**

 - Consolidate V&V outputs into a single ``docs/vnv/`` folder that each
   model change overwrites in place, so pull requests show before/after
   image diffs; ``vnv_growth_gif`` and ``vnv_compare`` now default to
   writing there

**v0.1.0 - 08/29/26**

 - Port simulation to the modern vivarium suite (vivarium 4.x, pandas 3.x)

   - Register particle columns via ``builder.population.register_initializer``
   - Route all particle state writes through ``Particle3D.update_particles``
     (vivarium 4 restricts column writes to the component that created them)
   - Use explicit attribute lists for all population view reads
   - Replace removed ``requires_columns``/``requires_values`` arguments with
     ``required_resources``
   - Fix force-cache index alignment so forces are attributed to the correct
     particles once some particles are frozen
   - Fix ``PathFreezer.get_population`` to look up KDTree neighbors positionally
   - Skip the end-of-simulation visualizer event loop in headless runs
   - Require Python 3.11+ (pandas 3); move unused GBD artifact tooling
     dependencies to the ``data`` extra

**v0.0.0 - mm/dd/yy**

 - Initial release

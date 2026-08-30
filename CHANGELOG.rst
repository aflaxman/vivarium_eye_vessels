**v0.4.0 - 08/30/26**

 - Paired arterial and venous trees (realism roadmap idea 3)

   - New ``vessel_type`` particle column: alternating artery/vein roots at
     the seed circle (arteries at ``particles.artery_caliber_ratio`` of the
     vein caliber), inherited through freezing, splitting, and DLA
   - ``FrozenRepulsion`` scales artery-vein repulsion by
     ``cross_type_factor`` so the trees interdigitate
   - ``PerfusionDemand`` is type-aware: tissue needs arterial supply and
     venous drainage separately, keeping the trees in balance
   - ``PathSplitter`` re-sprouts any tree whose active tips all died,
     instead of only when every tip in the simulation was gone
   - Visualizer and growth GIF color arteries red and veins blue
   - V&V: per-tree coverage fractions and the arcade A:V caliber ratio in
     metrics.json; bifurcation angles and junction exponents now measure
     only true bifurcations (continuation children excluded)

**v0.3.0 - 08/30/26**

 - Hypoxia-driven growth via space colonization (realism roadmap idea 2)

   - New ``PerfusionDemand`` force component: a lattice of tissue demand
     sites inside the ellipsoid; hypoxic sites (no frozen vessel within
     ``perfusion_radius``) each recruit their nearest growth tip, pulling
     tips toward unperfused territory until the tissue is colonized
   - ``model_spec.yaml`` enables it (``magnitude`` 0.3), deepens branching
     (``max_depth`` 4), and raises ``path_extinction.force_threshold`` to
     1.2 so tips survive long enough to colonize distant tissue
   - New ``perfused_fraction`` V&V metric; at 800 steps coverage rose from
     68.8% to 97.5%, vessel area density from 3.93% to 7.81%, and skeleton
     fractal dimension from 1.21 to 1.41 (HRF reference: 1.35 +/- 0.02)
   - Replace deprecated ``skimage.morphology.binary_dilation`` usage

**v0.2.0 - 08/29/26**

 - Vessel calibers via Murray's law (realism roadmap idea 1)

   - New ``radius`` particle column: roots seeded from
     ``particles.root_radius``, continuations tapered by
     ``path_freezer.radius_taper``, DLA attachments given
     ``path_dla.attach_radius``
   - ``PathSplitter`` assigns daughter calibers by Murray's law and derives
     branch angles from the radii via the minimum-work principle
     (``murray_exponent``, ``flow_asymmetry``, ``min_radius``)
   - Visualizer, growth GIF, and V&V rasterization draw true calibers
   - New V&V metrics: vessel area density and fitted junction exponents;
     caliber diagnostics overwrite ``docs/vnv/`` in place

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

**v0.13.0 - 08/31/26**

 - Comb-like side branching off the arcades (monopodial branching), with
   a branch-spacing calibration target

   - Visual review showed the branches off the arcades far too sparse:
     real arcades are monopodial -- a trunk that keeps nearly its own
     caliber and sheds small side branches at short, regular intervals
     -- while the splitter only did near-symmetric dichotomous forks,
     and the caliber cadence made wide tips branch *rarely*
   - New metric ``wide_junction_spacing_px``: skeleton distance between
     branch points along wide (>4 px) vessels, counting connected
     junction clusters once; measured identically on simulated rasters
     and HRF masks. HRF: a branch point every 22.7 px of wide skeleton
     (sd 1.77); the previous model: every 38.6 px. New calibration
     target with the across-mask sd as scale
   - New ``PathSplitter`` comb mode (``side_branch_flow``,
     ``side_branch_radius``, ``side_branch_probability``; off by
     default): parents wider than the reference split at their own
     emission probability -- exempt from the cadence damping -- and
     asymmetrically: the trunk continues at ~96% of its caliber while
     the tooth takes the Murray caliber for a ~10% flow fraction,
     leaving near-perpendicularly on a random side (both consequences
     of the minimum-work bifurcation relations already in the code)
   - The comb initially collapsed the network: teeth spawn inside the
     frozen-repulsion field of their own trunk and neighbours (radius
     0.15 units = 38 px, wider than the comb spacing itself) and the
     stacked forces drove every tip extinct. Real vessels pack a branch
     every ~23 px, so ``frozen_repulsion.interaction_radius`` drops to
     0.12; the extinction threshold is untouched
   - Fit: emission 0.65 per split round gives spacing 26.0 px (HRF 22.7)
     and calibration total 47.3, vs 135.7 for the previous spec under
     the corrected objective (the spacing term alone was 80.7); 0.8
     packs tighter (19.9 px, 50.5) but collapses a held-out seed, while
     0.65 holds on both held-out seeds
   - Three new tests: cadence exemption for side-branching trunks,
     asymmetric near-perpendicular tooth geometry, and the spacing
     metric on a synthetic comb

**v0.12.0 - 08/31/26**

 - Straighten the large vessels (caliber-dependent steering stiffness,
   an adaptation growth cap, and a fundus-faithful comparison)

   - Visual review showed the wide vessels meandering and curling at the
     zoomed-out scale, unlike the smooth HRF arcades. Two mechanisms were
     responsible: caliber-blind steering (arcade tips took the same
     random kicks as capillary tips, so their heading random-walked) and
     an anastomosis-shortcut runaway in the flow remodeler (an
     artery-vein capillary bridge sees enormous shear, thickens, gains
     conductance as r^4, draws more flow, and runs away to arcade
     caliber -- promoting curly capillary paths into fat vessels; on the
     fit seed only 376 of 2,516 wide particles were true arcades, and
     191 sat in the capillary-only deep plexuses)
   - New ``particles.noise_caliber_reference/_exponent``: tips wider
     than the reference caliber have their random steering attenuated by
     ``(reference/radius)^exponent`` (spec: 0.004 and 1.5), so arcades
     hold their heading while capillaries wander; exponent 0 is the
     legacy caliber-blind steering, and the OU disease dial composes
     with the attenuation unchanged
   - New ``flow_remodeler.max_adapted_radius`` (spec 0.006): shear-driven
     thickening saturates at venule caliber; segments born wider (Murray
     splits off the arcades) are untouched and can only taper. The
     default equals ``max_radius``, preserving the previous behavior
   - The HRF comparison now rasterizes the superficial (layer 0)
     projection only: fundus photographs do not see the deep
     capillary-only plexuses (OCTA does -- ``docs/vnv/plexus.png``), so
     scoring them against fundus masks inflated the capillary share.
     Tree- and graph-based metrics still use the full 3D network
   - ``initial_circle.n_vessels`` 4 -> 6: with the runaway capped, the
     missing wide-vessel mass was real arcade mass; three artery/vein
     pairs at the disc close most of the superficial density gap
   - New ``wide_tortuosity_q90`` calibration target (HRF 1.11, one-sided)
     guards against meandering wide vessels, though skeleton-branch
     tortuosity only weakly sees the long-wavelength curl (branches are
     chopped at junctions); the steering exponent and adaptation cap are
     therefore pinned during refits rather than left to the score
   - Two new tests: adaptation growth never exceeds the cap;
     a caliber-stiffened run lays straighter arcades than exponent 0

**v0.11.0 - 08/31/26**

 - Surface the diameter composition and fit it toward the HRF fractions

   - The comparison figure gains a fourth row: a per-branch diameter
     histogram (sim vs HRF, log-x) and grouped composition bars for the
     <=2 px / 2-4 px / >4 px strata, by branch count and by skeleton
     length, with the shares printed on the stratified panels;
     ``metrics.json`` records the diameter summaries and per-stratum
     length shares
   - The wide (>4 px) branch share joined the calibration targets;
     together with the capillary share this pins the whole three-bin
     composition (HRF: 6/59/35% capillary/mid/wide)
   - New ``flow_remodeler.adaptation_deadband``: segments whose shear is
     within this factor of their tree's median don't adapt (default 1.0
     = previous behavior). Motivated by the per-target calibration logs:
     shear adaptation polarizes calibers away from the median --
     below-median twigs thin to the floor while above-median vessels
     thicken -- hollowing out the mid-caliber stratum that dominates
     real networks
   - Composition fit: ``adaptation_rate`` 0.10 -> 0.05,
     ``shear_threshold_fraction`` 0.5 -> 0.65,
     ``capture_radius`` 0.045 -> 0.035 moves the composition
     19/15/66% -> 21/21/58% and the score 40.2 -> 34.8 (scores now
     include the wide-share target, so they are not comparable to
     v0.10.0's 24.7); on the fit seed reverting any single knob of the
     trio is worse, and the trio also wins on a held-out healthy seed
     (54.8 -> 40.5) but amplifies the pre-existing degeneracy of seed 42
     (165 -> 326), keeping the multi-seed objective on the roadmap
   - The residual composition gap to HRF is structural -- Murray
     bifurcations plus geometric taper transit the mid-caliber band in a
     couple of generations -- and is called out as mechanism work in the
     roadmap

**v0.10.0 - 08/31/26**

 - Calibrate the healthy model against real-data targets (realism roadmap
   idea 8, healthy fit)

   - New ``vnv_calibrate`` CLI: a formal objective over eight validation
     targets (HRF-derived skeleton/area density, fractal dimension, and
     tortuosity; KS distance to the HRF segment-length distribution;
     capillary branch share; clinical A:V ratio; one-sided full
     perfusion), each a squared z-like deviation, fitted by coordinate
     descent over the eight most influential knobs with every 800-step
     evaluation logged
   - ``metrics.json`` now records the calibration score of every standard
     run (``calibration`` block + figure headline), tracked across model
     versions like every other metric
   - The first fit (22 evaluations) moved three knobs
     (``caliber_cadence_exponent`` 0.6 -> 0.75, ``adaptation_rate``
     0.15 -> 0.10, ``perfusion_demand.magnitude`` 0.3 -> 0.35) and cut
     the score 38.7 -> 24.7: KS to the HRF length distribution
     0.171 -> 0.074, capillary share 22.5% -> 19.1%, fractal dimension
     1.43 -> 1.41, area density 14.2% -> 13.6%, A:V ratio and tortuosity
     holding; skeleton density traded back 2.69% -> 2.43%

**v0.9.0 - 08/30/26**

 - Smooth, controllable tortuosity via Ornstein-Uhlenbeck steering
   (realism roadmap idea 7)

   - New ``wx/wy/wz`` particle columns: the autocorrelated random
     component of each tip's steering, updated as an AR(1) with
     correlation time ``particles.noise_persistence_time`` (days) and the
     same stationary spread as the legacy uniform kick, so the process
     degenerates exactly to white noise at one-step persistence and
     ``noise_persistence_time: 0`` preserves the legacy behavior
     bit-for-bit
   - The healthy spec keeps ``noise_persistence_time: 0``: the sweep
     showed healthy tortuosity was already at the HRF target (branch
     median 1.004 vs 1.000), and any persistence raises tortuosity
     (1.011 -> 1.049 path median at 0.25 days) as a package with
     capillary dropout and hypoperfusion -- a DR-like disease dial, not
     a healthy-network knob

**v0.8.0 - 08/30/26**

 - Layered plexuses (realism roadmap idea 6)

   - New ``layer_id`` particle column: each vessel's home plexus,
     inherited through freezing, splitting, and DLA
   - New ``PlexusLayers`` component: a damped Hookean spring holds each
     active tip near its layer's z-plane (undamped, tips oscillate,
     saturate the terminal-velocity clamp, and go extinct), and
     capillary-caliber tips dive one layer deeper with
     ``dive_probability`` — the frozen transit trail is the diving
     vessel; wide vessels never dive, so deeper plexuses are
     capillary-only
   - Force rebalance for planar crowding: ``path_extinction``
     ``force_threshold`` 1.2 -> 2.0, anastomosis ``capture_radius``
     0.03 -> 0.045 (reaches across a plexus gap)
   - V&V: new ``docs/vnv/plexus.png`` with OCTA-style en-face slabs per
     layer plus an x-z cross-section; ``plexus_layers`` block in
     ``metrics.json`` (per-layer counts, calibers, z-adherence, diving
     vessels)
   - At the standard 800-step run: plexuses hold 5,079 / 931 / 1,624
     segments with 173 diving vessels; headline metrics intact (97.9%
     perfused, skeleton density 2.69%, area density 14.2%, fractal
     dimension 1.43, arcade A:V ratio 0.68) and anastomosis loops rose
     from 21 to 39

**v0.7.0 - 08/30/26**

 - Flow-based remodeling and pruning (realism roadmap idea 5)

   - New ``FlowRemodeler`` component: solves Poiseuille flow on the frozen
     graph as a resistor network (conductance r^4/L per segment and
     anastomosis bridge, artery/vein roots as fixed-pressure terminals, a
     per-node tissue leak) via one sparse Kirchhoff solve per
     ``remodel_interval``
   - Wall shear |Q|/r^3 drives remodeling: low-shear *terminal* segments
     (degree-1 graph ends, so the network is never cut mid-branch) are
     pruned and recycled into the free particle pool, stamping the
     previously unused ``unfreeze_time`` column; every caliber drifts by
     ``adaptation_rate`` toward the median-shear radius (trunks thicken,
     twigs thin)
   - New ``docs/poiseuille_flow.md``: the physics from first principles —
     the Hagen-Poiseuille law, the r^4 punchline, the Ohm/Kirchhoff
     circuit analogy, wall shear as the vessel's own sensor, and how
     Murray's law emerges from uniform shear
   - V&V: ``n_pruned`` and ``wall_shear`` summaries in ``metrics.json``,
     pruned count in the figure headline
   - Shear targets are per tree (arteries genuinely run at higher wall
     shear than veins) and the depth-0 arcades are exempt from adaptation
     (their calibers are upstream boundary conditions); both were needed
     to keep the clinical A:V caliber ratio, which holds at 0.66
   - At the standard 800-step run: 3,281 segments pruned, capillary-scale
     branch share 66% -> 27% (HRF 6.5%), vessel area density 7.2% -> 15.1%
     (HRF 11.9%), fractal dimension 1.49 -> 1.43 (HRF 1.35), tissue 98.9%
     perfused; the junction exponent becomes emergent (median k ~ 2.2,
     close to measured human values) instead of the imposed 3.0, and the
     aggregate length median drifts to 29 px (HRF 22) as clutter is
     removed

**v0.6.0 - 08/30/26**

 - Anastomosis: capillary loops between the trees (realism roadmap idea 4)

   - New ``PathAnastomosis`` component: an active capillary-caliber tip
     that comes within ``capture_radius`` of the *other* tree's
     capillary-caliber frozen segments fuses onto it (with a per-step
     ``probability``), freezing the tip and recording the join in a new
     ``anastomosis_id`` particle column — capillaries join capillaries,
     tips don't fuse into trunks
   - ``tree_edges``, the rasterizer, the visualizer, and the growth GIF
     draw the bridges (violet), so the network is rendered and measured
     as the perfusable graph it now is
   - New ``graph_cycles`` V&V metric (E - N + C over parent and
     anastomosis edges) counts independent loops, plus ``n_anastomoses``
     in ``metrics.json``; both were structurally zero for the tree-only
     model
   - Enables flow-based remodeling (roadmap idea 5): Poiseuille flow now
     has closed circuits to solve on
   - ``metrics.json`` now records wall-clock runtime (setup, simulation,
     steps/second) so speed regressions are tracked alongside the network
     metrics

**v0.5.0 - 08/30/26**

 - Caliber-dependent branching cadence (realism roadmap idea 1b)

   - ``PathSplitter`` scales each tip's split probability by
     ``(min_radius / radius) ** caliber_cadence_exponent``, so wide trunks
     run long between branch points while capillary-caliber twigs branch at
     the full ``split_probability`` (exponent 0 restores the old
     caliber-independent cadence)
   - ``model_spec.yaml`` pairs exponent 0.6 with a faster base cadence
     (``split_interval`` 15, was 30) so twigs branch more densely while
     trunks branch less
   - Targets the segment-length distribution — the clearest visual gap
     against the HRF reference: at 800 steps the KS distance to the HRF
     log-length distribution fell from 0.099 to 0.053, the median matches
     HRF exactly (22 px), the shortest-length spike dropped from 2.0x to
     1.6x the HRF density, and skeleton density rose from 3.03% to 3.39%
     (HRF: 3.21%)
   - New tree-based ``tree_segment_length`` metric in ``metrics.json``:
     inter-branch-point distances in simulation units, independent of
     rasterization
   - New diameter-stratified segment lengths: per-branch diameter is
     recovered from the distance transform along the skeleton (so real
     masks without explicit calibers are measured the same way as the
     sim), and lengths are compared within capillary / mid / wide strata —
     a new figure row and ``branch_length_by_diameter`` in ``metrics.json``
   - Side benefit: the arcade A:V caliber ratio moved from 0.84 to 0.67,
     matching the clinical target

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

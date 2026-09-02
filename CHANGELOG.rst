**v0.21.0 - 09/02/26**

 - Hypoxic survival relief: an extinction-threshold pipeline (negative
   result, default-off)

   - ``PathExtinction`` now reads its per-tip threshold from a
     ``particle.extinction_threshold`` value pipeline (base
     ``force_threshold``), and ``PerfusionDemand`` registers a modifier
     that raises it by ``perfusion_demand.survival_factor`` for tips
     within ``perfusion_radius`` of a hypoxic site of their own type --
     hypoxia-driven VEGF is a survival signal for tip cells as well as a
     chemoattractant. This is the extinction-side lever the eleventh
     pass pointed at: a sprout born beside a dense frozen frontier gets
     a new ``path_id``, so the frontier's repulsion hits it at full
     strength from step one and pushes it over the threshold before it
     can escape into the tissue that recruited it. The pipeline wiring
     avoids an import cycle and keeps the mechanism with the demand
     field it depends on. 1.0 is bit-for-bit legacy; two unit tests
   - It ships OFF. Sweeping 1.5/2/3 on the calibration seeds: the fit
     seed is bit-identical across all three (51.0 -> 61.8, still fully
     perfused -- the relief is effectively binary, the same tips are
     saved at every factor), seed 7 improves enormously (168.5 -> 70.1 /
     67.0 / 67.0 at full perfusion) and seed 909 collapses at every
     factor (125.4 -> 967.2 / 884.2 / 884.2, perfusion 0.99 -> 0.44 /
     0.47; 3-seed means 366.3 / 337.7 / 337.7 against 115.0). Saved tips
     on 909 do not build usable vessel: the network shrinks 11098 ->
     3240-4035 particles
   - Third lever in a row -- after the ninth pass's growth gates and
     the eleventh's prune grace -- on a third subsystem to show the same
     seed-7-vs-909 trade. The roadmap's twelfth pass draws the
     conclusion: post-hoc per-tip relief cannot fix a trajectory decided
     early, and the remaining seed variance most plausibly sits in the
     random initial root headings; a radially-outward disc emergence (as
     real arcades do) is the next candidate

**v0.20.0 - 09/02/26**

 - Prune grace period (a negative result, default-off) and a DRY pass
   over the stacked increments

   - ``flow_remodeler.prune_grace_days``: terminal segments frozen less
     than this long ago are not pruned, so a young sprout gets time to
     connect and earn its shear before regression judges it. Built to
     close the tenth pass's remaining gap (two held-out seeds settle at
     ~0.88 perfusion because pruning consumes what the held wave front
     builds). Age-based sparing was chosen over a position band around
     the front: one knob, no cross-component coupling, and it also
     protects the wave's behind-front re-sprouts. Unit-tested
   - It ships OFF. Sweeping 2.5/5/10 days on the calibration seeds: 2.5
     is worse on every seed (fit 51.0 -> 60.1, seed 7 168.5 -> 241.7,
     seed 909 125.4 -> 154.8); 5 and 10 rescue seed 7 spectacularly
     (-> 80.3 / 86.0 at full perfusion, 22-26k particles) but halve the
     fit seed's network (-> 121.5 / 126.6, perfusion 0.87-0.89,
     ``capillary_share`` 0.10 -> 0.35) and collapse seed 909 (125.4 ->
     580.0 / 641.4, perfusion 0.59 / 0.55); 3-seed means 152.2 / 260.6 /
     284.7 against the incumbent 115.0. Mechanism: pruning declutters
     as well as remodels -- a spared dead-end stub stays frozen and acts
     as a repulsion source that chokes the tips around it; whether the
     stubs become vessels or clutter is seed-specific. The same per-seed
     reshuffle signature as the ninth pass's local rules, so the wave's
     0.87-1.00 band stands as the shipping state
   - Behavior-preserving consolidation, verified bit-for-bit over 150
     steps: ``PathSplitter.eligible``/``commit`` replace three copies of
     the depth-ceiling + crowding-gate + write sequence;
     ``PerfusionDemand.vessel_distances`` replaces the nearest-vessel
     KDTree query duplicated in ``hypoxic_sites`` and the wave's
     ``served_fraction`` (the wave drops its own freezer handle);
     ``DevelopmentalWave.disc_distance`` replaces two inline norms;
     ``vnv.simulation.with_seed``/``build_from_spec`` replace the
     seed-override and write-then-build boilerplate in ``vnv_calibrate``
     and ``vnv_contact_sheet``. ``docs/vnv`` is unchanged because the
     spec's dynamics are unchanged

**v0.19.0 - 09/02/26**

 - Developmental wave: a closed-loop growth front, on by default

   - New ``DevelopmentalWave`` component: vascularization advances as a
     wave from the optic disc, as in real retinal development.
     ``PerfusionDemand`` exposes only the demand sites within
     ``radius + lookahead`` of the disc, so every tip chases the same
     expanding ring; the front advances by ``wave_speed`` per step only
     while the tissue behind it is served (``advance_threshold``); and a
     front held ``hold_resprout_steps`` re-sprouts the stalled tree from
     the frozen vessels nearest its unserved sites -- growth pressure
     applied exactly where growth is missing, which is what the ninth
     pass's untargeted tip floor lacked. The control variable is the
     goal itself (perfusion behind the front), not a local proxy
   - First mechanism in the project to improve every calibration seed
     at once: 3-seed mean 581.4 -> 115.0 (fit seed 47.0 -> 51.0 at full
     perfusion; seed 7 375.8 -> 168.5; chronic collapser seed 909
     1321.4 -> 125.4 at 0.99 perfusion). A ten-config sweep set the
     values: ``advance_rule: combined`` (per-type advance brakes healthy
     seeds on the artery tree's ~0.85 mid-field service ceiling and
     pruning grinds the held network), threshold 0.85 (0.75 outruns
     weak seeds, 0.95 makes strong seeds wait), and no wide-tip
     exemption (it breaks the front discipline)
   - Held-out contact sheet: the catastrophic collapse mode is gone --
     min perfusion across the four held-out seeds rises 0.40 -> 0.87,
     mean 0.83 -> 0.94, mean score 574.7 -> 255.3. The honest trade:
     strict >=95% reliability stays 2/4 because two previously-fast
     seeds settle at a stable ~0.87-0.88-perfused equilibrium (verified
     out to step 1200: the wave-front demand and shear pruning balance
     just below the bar) -- a complete connected network rather than
     legacy's collapse, but not full perfusion. The wave trades
     legacy's bimodal outcomes (0.99 or 0.40) for a compressed 0.87-1.00
     band; decoupling pruning from held-front regions is the natural
     follow-up
   - ``vnv_compare``/``vnv_contact_sheet`` artifacts regenerated under
     the new default; ``resprout_at`` on ``PathSplitter`` gives the
     wave its targeted re-sprout entry point; four new unit tests
     (74 total)

**v0.18.0 - 09/02/26**

 - Growth reliability: mechanisms, a held-out-seed V&V artifact, and an
   honest negative result

   - A field study on four fresh seeds (11/202/909/4242) found two of
     four networks stall mid-growth -- the collapse the roadmap's
     percolation analysis predicted. Three candidate fixes are now
     implemented and unit-tested; all three ship default-off, because
     every configuration that improved the calibration seeds hurt the
     held-out seeds (details below and in the roadmap's ninth pass)
   - Crowding gate (``path_splitter.max_crowding`` within
     ``crowding_radius``): skip a dichotomous split when the tip already
     has that many frozen neighbors, keeping branching at the growth
     front where daughters can survive; comb teeth are exempt (a
     trunk's own frozen trail dominates its neighbor count)
   - Self-healing growth front (``path_splitter.min_active_tips``,
     gated by new ``path_splitter.resprout_established_size``): a tree
     whose active-tip count thins below the floor re-sprouts from its
     own frozen vessels, but only once it is established -- early on,
     few tips is the natural state and topping up just crowds the disc
   - Balanced arterial inflow (``flow_remodeler.balanced_arterial_inflow``):
     split the total arterial inflow equally across artery roots as
     current sources instead of fixing every root pressure, so no
     arcade can starve its siblings (the rich-get-richer failure:
     the widest arcade takes the flow, shear adaptation amplifies it)
   - The negative result, on record so it is not re-litigated: the best
     multi-seed sweep config (gate 30 + tip floor 8, established at
     200) improved the 3-seed calibration mean 581.4 -> 431.8, but
     held-out reliability fell from 2/4 seeds perfused to 0/4.
     Attribution probes: the tip floor alone drops seed 11 from 0.996
     to 0.76 perfused and seed 4242 from 0.928 to 0.54; the gate at 30
     drops seed 202 from 0.987 to 0.79 (healthy interiors run 30-44
     frozen neighbors, so any gate low enough to fire also caps normal
     filling -- a gate of 60 never fires at all, bit-for-bit). The
     growth front is near-critical: these knobs reshuffle which seeds
     fill rather than making every seed viable
   - New ``vnv_contact_sheet`` CLI renders held-out-seed simulations
     (11/202/909/4242, none used in calibration) beside HRF expert
     masks, stamping each panel with perfusion and skeleton density and
     a [STALLED] flag; ``docs/vnv/contact_sheet.json`` records the
     seed-reliability fraction so reliability is tracked across
     versions like every other metric -- this artifact is what caught
     the non-generalization above

**v0.17.0 - 09/01/26**

 - Sub-teeth: extend the comb one caliber class down, with the
   robustness trade on record

   - Visual review still read the mid vessels as too long and meandering
     between branch points. Diagnosis: comb teeth ran ~90 px between
     splits (cadence damping) where HRF mid vessels branch every
     30-50 px, and artery teeth sat just below ``side_branch_radius``
     so they never combed. ``side_branch_radius`` 0.008 -> 0.006 lets
     teeth grow sub-teeth -- the hierarchical texture real mid-vessels
     show -- and the fit-seed network now pairs it with near-target
     skeleton density (2.9%) and full perfusion
   - The robustness trade is deliberate and documented: the 3-seed mean
     regresses (67.4 -> 165) because one seed in three stalls mid-growth
     under the denser branching. Nine stabilizer configs were swept
     multi-seed (extinction thresholds, repulsion radii, developmental
     gating, self-healing growth fronts, bigger particle pools -- the
     pool hypothesis was falsified bit-for-bit); none kept the sub-teeth
     look while protecting every seed. The growth front is
     percolation-like, and the roadmap names this the central open
     problem for the healthy model
   - Three new mechanism knobs from the investigation, all
     default-legacy and unit-tested: ``path_anastomosis.min_layer``
     (restrict capillary fusion to deeper plexuses),
     ``path_splitter.min_active_tips`` (a tree whose growth front thins
     re-sprouts to top itself up; 1 = legacy bootstrap-only), and
     ``path_splitter.side_branch_start_time`` (combs form after the
     arcades establish, as in retinal development)

**v0.16.0 - 08/31/26**

 - Match the length-weighted caliber profile ("length x width")

   - New metric ``skeleton_pixel_diameters``: the local (2 x EDT)
     diameter at every skeleton pixel, so its distribution is skeleton
     length by width with no binning into strata. The comparison figure
     replaces the branch-diameter histogram with this caliber profile,
     and ``ks_caliber_profile`` (KS between the sim's superficial
     profile and the pooled HRF profile) joins the calibration targets
   - The profile located the residual mismatch precisely: an excess at
     ~2.8 px -- vessels piling up exactly at ``max_adapted_radius`` --
     and a 3x deficit in the 5.5-7 px band real fundi carry. Vein comb
     teeth are born at 4.8 px wanting to thicken (their shear sits above
     the tree median), but the cap only lets them shrink: 76% of depth-1
     particles had been ground to <= 3.2 px, median exactly at the cap
   - Fit (multi-seed): ``adaptation_deadband`` 1.0 -> 2.0 (moderate-
     shear segments keep their born caliber instead of being ground
     toward the cap) plus ``side_branch_flow`` 0.1 -> 0.15 (teeth carry
     15% of trunk flow, so vein teeth are born at ~5.4 px, inside the
     missing band). 3-seed mean 79.7 -> 67.4, caliber-profile KS
     improves on every seed, obtuse-angle share halves. Raising the cap
     itself (0.006 -> 0.010, with or without the deadband) re-triggers
     the shortcut runaway (means 244 / 161) and stays rejected.
     Deadband-alone edges the aggregate mean (64.2) by winning back
     seed-42 density terms, but loses the caliber profile and angle
     terms this pass targets; the trade is recorded here deliberately

**v0.15.0 - 08/31/26**

 - Unimodal bifurcation angles: measure fundus-visible junction geometry
   and pin it with targets

   - The bifurcation-angle histogram was bimodal (a peak near 75 degrees
     plus a second mode at 100-150). Decomposing junctions by provenance
     showed the second mode is almost entirely deep-plexus capillary
     junctions (layer-2 share of obtuse angles: 0.80; superficial: 0.08)
     -- polygonal-mesh T-junctions that fundus photographs cannot see,
     and that real deep plexuses genuinely have. The visible (arcade and
     comb) junctions were already unimodal at ~76 degrees with zero
     obtuse share
   - The comparison figure's angle panel and Murray-exponent inset now
     measure the superficial tree only, consistent with the raster
     (``metrics.json`` keeps the all-layer summary alongside)
   - Two new literature targets on the superficial tree:
     ``bifurcation_angle_median`` (77 +/- 5 degrees) and one-sided
     ``bifurcation_obtuse_share`` (>100 degrees; target 0.05). The
     current spec sits at median 75-79 and obtuse 0.09-0.14 across the
     three calibration seeds
   - New ``perfusion_demand.caliber_reference/_exponent`` knob
     (default off, tested): attenuates hypoxia chemotaxis on wide tips.
     Probed as a mechanism for the residual visible splay and REJECTED
     by the multi-seed objective -- the pull on wide tips is what drives
     the arcades outward, and attenuating it stalls colonization
     (3-seed mean 72.3 -> 2563 at exponent 1.0). The knob stays for
     disease-phenotype work

**v0.14.0 - 08/31/26**

 - Multi-seed calibration objective, and a multi-seed refit

   - ``vnv_calibrate`` gains ``--seeds A,B,C``: the objective becomes
     the MEAN score across the listed seeds (one simulation per seed
     per evaluation, per-seed breakdowns logged), so a config that
     collapses on any seed loses to one that is merely mediocre
     everywhere -- the robustness the single-seed fits kept missing
   - The comb-era candidates were re-judged on seeds 123456/7/42:
     every config has exactly one weak seed, and the single-seed
     winner's weak seed was hiding real fragility. The best 3-seed
     mean is ``plexus_layers.dive_probability`` 0.05 -> 0.04
     (mean 114.1 -> 70.6; per-seed 47/74/221 -> 42/121/49), which
     also improves the fit seed itself
   - Two new tests: the multi-seed objective is the per-component mean,
     and a config that collapses on one seed loses to a steady one

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

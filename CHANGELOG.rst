**v0.30.0 - 09/06/26**

 - The capillary bed: capillary-scale hypoxia, sprouting from vessel walls

   - **Why.** The model's finest vessels were 18 um wide and about 0.7 mm
     apart (the arteriole-scale demand lattice declares tissue within
     0.68 mm of a vessel perfused); a real superficial plexus has 8 um
     capillaries 60-80 um apart. That one coarse scale sat behind the
     voids in the contact sheet, the empty perifovea of the previous
     round, and the deep plexus's coarseness. The ROSE-1 labels measure
     it: intervessel distance 79 +/- 20 um and skeleton density 9.3 +/-
     2.4 mm/mm2 in the 3 x 3 mm macular window (outside the FAZ),
     against the model's 330-400 um and 2.1-2.3 mm/mm2
   - **Measurement.** ``metrics.capillary_statistics`` (twice the mean
     distance from a non-vessel pixel to the nearest skeleton pixel, and
     skeleton length per area, both outside the FAZ) on the ROSE-1 SVC
     expert labels (``reference_data.fetch_rose_labels``) and on the
     model's OCTA window, which now draws every segment at least a pixel
     wide because OCTA images flow, not caliber. Two new targets,
     ``octa_intervessel_um`` and ``octa_skeleton_mm_per_mm2``, derived by
     ``--derive-targets``. The labels omit some perifoveal capillaries, so
     the reference spacing reads a little wide and the density a little
     low; the model is measured without that omission
   - **Model.** ``CapillaryBed``: inside a disk around the fovea (2.7 mm,
     the OCTA window with margin), for each plexus layer, a 90 um lattice
     of tissue sites; a site with no frozen vessel of its layer within
     90 um is hypoxic. Every 5 steps up to 40 hypoxic sites behind the
     developmental front recruit a sprout from the nearest vessel wall of
     their layer within 225 um: an 8 um tip aimed at the site (angiogenic
     sprouting from a wall, not tip splitting). Capillary tips are pulled
     toward their layer's hypoxic sites, tolerate three times the
     extinction force, and freeze where no hypoxic site remains within
     180 um. Loops close through PathAnastomosis; dead-end sprouts carry
     no flow and are pruned. Hooks so the rest of the model treats them
     as capillaries: they never split (a Murray split would hand them
     wider daughters), never dive to another layer, feel FrozenRepulsion
     only within 90 um instead of 540 um, and neither set a tree's median
     shear nor adapt their caliber in the FlowRemodeler. Capillaries are
     told from arterioles by caliber alone, so they sit strictly below the
     remodeler's adaptation floor (a cutoff at the arteriole floor, 0.002,
     also dropped the thinned arterioles and shrank the baseline network
     by a third). An 8 um vessel is
     below half a pixel at fundus scale, so the HRF comparison does not
     see them except through what they do to the arterioles
   - **Coupling, learned the hard way.** With the bed in, the macula's
     arterioles thinned to a third (spec seed: 1361 -> 374 arteriole-scale
     particles within 2.7 mm of the fovea) and the fundus-visible macular
     clear radius grew to 1.0 mm, because the rest of the model treated
     capillaries as vessels like any other. Three more hooks say what a
     capillary is not: it does not perfuse a PerfusionDemand site (a
     capillary carries blood only when an arteriole feeds it, so the coarse
     lattice keeps recruiting arterioles into tissue the bed has reached;
     ``perfusion_demand.min_radius``), it does not fence out an arteriole
     tip (FrozenRepulsion ignores frozen capillaries for wider tips), and
     it does not count as crowding at a PathSplitter branch point, and a
     precapillary tip does not end in it (PathAnastomosis fuses only
     capillary tips onto capillaries). The bed is also kept out of the
     FlowRemodeler's solve (``flow_remodeler.capillary_radius``): confined
     to the macula it was a local sink that stole flow from every other
     branch -- arteriole pruning rose a quarter and on seed 7 an arcade
     curled around the macula and closed on itself, leaving a 1.75 mm
     void. The bed regresses its own dead ends instead (a stopped sprout
     that closed no loop unwinds after 2.5 days, tip first). Junction
     statistics compared against fundus references skip capillary sprouts,
     as a fundus does (``metrics.CAPILLARY_RADIUS_UNITS``). Two more couplings surfaced on the seed set: the arteriole
     trees' re-sprouting (the tip floor in PathSplitter and the wave's
     stall response) picked frozen capillaries as walls to sprout from,
     and the bed's churn through the shared particle pool made the
     splitter skip split rounds -- capillary walls are now excluded and
     the bed tops the pool up ahead of both consumers
   - **Results.** The bed meets the capillary-scale targets on every
     seed: intervessel distance 66-71 um (ROSE 79 +/- 20) and skeleton
     density 9.7-10.7 mm/mm2 (9.3 +/- 2.4) in the OCTA window, with a
     capillary ring on the FAZ edge (clear radius 0.14-0.31 mm, target
     0.29 +/- 0.07; the foveal exclusion radius returns to 0.065 now that
     capillaries, not arterioles, define the zone). Eight seeds
     (7/42/909/2024/123456/31/77/5150): mean 53.9 (97/50/25/28/96/34/58/
     45) under an objective with two new terms the previous round could
     not score. Every term but one sums to 29.3 against 27.9 before; the
     one is the fundus-scale macular clear radius, 22 points of the mean
     on its own. It is the largest vessel-free disk within 1.5 mm of the
     window center, a maximum statistic that a single displaced arteriole
     moves by half a millimeter: the same seed reads 0.64 mm under one
     rule for withdrawing capillary tips from the FAZ and 1.30 mm under
     the other, with the macula's own arterioles unchanged (357 within
     2.7 mm of the fovea against 346 before the bed). Its eight-seed mean
     is 0.59 mm for v0.29, 0.69 for a bed-off control, 0.73 and 0.85 for
     the two bed runs (HRF 0.53 +/- 0.09), so the bed may open the
     macula's visible voids a little; the seed set cannot say. Paired
     perfusion dips below 0.9 on two seeds (42: 0.885, 5150: 0.879); the
     FAZ clear radius averages 0.23 mm (target 0.29). With the bed
     switched off every seed reproduces the previous round exactly. A
     control settles what the differences are worth: the bed off and only
     the particle pool changed from 500 to 501 -- which displaces the
     random draws the arteriole tree makes and nothing else -- scores
     45.0 on the terms v0.29 scored 29.9 on (per seed 33/50/52/27/56/52/
     62/28 against 44/18/42/19/25/26/29/36), paired perfusion 0.88-0.99,
     macular clear radius up to 0.91 mm. The eight-seed mean moves by 15
     points and single seeds by 30 under a change with no physics in it,
     so differences of that size between rounds are noise, and so was
     some of v0.29's fit. The seed set is large enough to find a 30-point
     collapse mode (the twentieth pass did), not a 6-point shift.
     A whole-field bed (``region_radius`` 0) runs about three hours a
     seed against 25 minutes for the macular one and is the natural next
     step once the flow solve can carry it
   - Held-out seeds 11/202/909/4242 all colonize 100% of the tissue
     (4/4) with arterial supply 1.00/0.99/0.93/1.00 (0.99/1.00/1.00/0.99
     in v0.29); their macular clear radii read 0.93/0.77/0.63/0.80 mm
     (1.07/0.56/0.62/0.98 before)
   - **Eccentricity, checked and deferred.** Whether the eye's growth
     earns a round was tested on HRF: skeleton density falls from 5.0%
     within 1.5 mm of the disc to 2.8% at 4 mm and the temporal side is
     1.5x denser than the nasal, while the model is flat (3-4% at every
     distance) and symmetric (ratio 0.9-1.1). Both gaps are as consistent
     with the missing perifoveal supply as with non-uniform growth, so
     they are recorded as targets for later rather than attributed

**v0.29.0 - 09/05/26**

 - The fovea: OCTA targets, and an exclusion that did not exclude

   - **Data.** The ROSE OCTA dataset (Retinal OCT-Angiography vessel
     SEgmentation; registration-gated, so ``reference_data.
     fetch_rose_images`` reads it from the cache directory rather than
     downloading it). ROSE-1 has 39 fovea-centered 3 x 3 mm scans of the
     superficial vascular complex (304 px, 9.9 um/px) with expert labels
     for large vessels and capillary centerlines, from a mix of healthy
     and Alzheimer's eyes. The labels omit many perifoveal capillaries,
     so the FAZ is read from the angiograms themselves
   - **Scale.** One simulation unit is 4.5 mm (``metrics.MM_PER_UNIT``):
     the model's disc-to-fovea distance (1.05 units) is the anatomical
     4.76 mm, the root vein caliber then reads 153 um and the HRF field
     15 mm, consistent to ~5%. The fundus raster is 17.6 um/px
   - **Measurement.** Two macular targets, both the radius of the largest
     vessel-free disk near the fovea (``metrics.clear_radius_px``),
     because an inscribed disk cannot leak through a gap the way a region
     can. At fundus scale, ``macular_clear_radius_mm`` on the skeleton
     within 1.5 mm of the image center: HRF 0.53 +/- 0.09 mm, the macula
     a photograph shows clear because its capillaries are too fine to
     see. At OCTA scale, ``faz_radius_mm`` on the vascular signal (the
     en-face image smoothed at 40 um, thresholded at half its mean --
     the mean, so a sparse network drawn on black still has a threshold)
     within 0.3 mm of the fovea: ROSE-1 0.29 +/- 0.07 mm. The simulation
     is read the same way, its superficial plexus rasterized on a 3 x 3
     mm window at the fovea (``metrics.octa_window``). The avascular
     region's area, the OCTA literature's FAZ area, is kept as a
     diagnostic: read from label masks it leaks into the intercapillary
     spaces (0.9 mm2 against 0.37 from the images, no closing radius
     repairs it), and read from the model's sparse window it is
     meaningless (4.2 mm2 with a trunk through the fovea)
   - **Finding.** The model had no fovea. The exclusion cylinder (radius
     0.4 units = 1.8 mm, spring 3, and a ``height`` key the component
     never read, so it spanned every plexus) was porous: 950 frozen
     particles sat inside it, the network's density there was only 13%
     below the surround, and an arcade trunk ran straight through the
     foveal center. The macula's largest visible void was 0.93 mm in
     radius (HRF 0.53) and lay away from the fovea; the FAZ radius read
     0.33 mm only because a gap happened to sit beside the trunk. Spec
     seed under the new terms: 68.6, of which the two macular terms 53.5
   - **Model.** The exclusion made real and caliber-aware. ``cylinder_
     exclusion`` gains ``wide_radius`` / ``wide_min_radius``: tips at or
     above 0.0025 (arterioles and venules) are held out to 0.10 units
     (0.45 mm), capillary tips to ``radius`` 0.05 (0.23 mm), with a
     spring of 30 so a trunk under the arcade guidance and hypoxia pulls
     stops within 0.02 units of the boundary (at 10 it leaks: FAZ radius
     0.13 mm). The dead ``height`` key is gone from the spec. The demand
     lattice no longer has sites inside the capillary radius: the fovea
     is fed by the choroid, so no tip is recruited into it. Spec seed
     68.6 -> 24.9 (FAZ clear radius 0.33 mm, macular clear radius 0.67
     mm); ``vnv_compare`` now writes ``docs/vnv/macula.png``
   - Eight seeds (7/42/909/2024/123456/31/77/5150): mean 29.9 (per seed
     44/18/42/19/25/26/29/36); the non-macular terms sum to 25.0 against
     25.6 in the previous round, so the fovea costs the rest of the fit
     nothing. FAZ clear radius 0.39 mm on average (0.30-0.54: on two
     seeds the capillary annulus stays empty and the zone reads the wide
     radius), macular clear radius 0.59 mm (0.44-0.67), paired perfusion
     0.96 (seed 7 falls to 0.91 with a trunk that used to feed the macula
     now routed around it)
   - Swept and rejected on the spec seed unless noted: radii 0.065/0.118
     (8-seed mean 45; the annulus between the two radii stays empty and
     the FAZ reads the wide radius, 0.45 mm on average); 0.05/0.118
     (47.5: seed 77 collapses to 140 with a 1.2 mm macular void); a
     single radius for every caliber (74 at 0.05, 132 at 0.065: with
     nothing to stop them short of the zone the arcades route around
     the whole macula and leave a 0.9-1.4 mm void); spring 10 (41.5,
     FAZ 0.13 mm: leaks); radii 0.04/0.09 (81); a denser demand lattice
     to fill the perifovea (``perfusion_radius`` 0.12: 43, 0.10: 127 and
     with ``site_spacing`` 0.05: 226 -- more demand than the front can
     serve, paired perfusion 0.85). The macular clear radius is a
     maximum statistic and swings 0.5-1.2 mm between near-identical
     variants on one seed; it needs the full seed set to rank anything
   - Held-out seeds 11/202/909/4242 all colonize 100% of the tissue
     (4/4) with arterial supply 0.99/1.00/1.00/0.99 (0.96/0.98/0.96/
     1.00 before); their macular clear radii read 1.07/0.56/0.62/0.98
     mm -- the network's own voids, on two of four seeds larger than
     the fovea
   - ``vnv_compare`` writes ``docs/vnv/macula.png`` (the macula at both
     scales, simulation beside HRF and ROSE); the contact sheet stamps
     each seed's macular clear radius. ``scripts/sweep/sweep_slurm.sh``
     runs a job file as a Slurm array

**v0.28.0 - 09/04/26**

 - Arcade geometry: the front's demand gate curled the trunks

   - **Measurement.** Three fundus-scale statistics of the wide (> 4 px)
     skeleton, read identically on HRF masks and on a fundus-sized window
     of the simulated raster (``metrics.arcade_geometry``), relative to
     the optic disc estimated from the image itself as the point the wide
     vessels' tangent lines converge on (``metrics.disc_center``; on the
     disc in all 15 HRF masks, within ~40 px of the true disc on the
     simulation): ``arcade_radial_alignment``, the mean |cos| between a
     wide vessel's tangent and the direction from the disc (HRF 0.81 +/-
     0.03); ``arcade_reach_px``, the mean distance of wide skeleton from
     the disc (247 +/- 14 px); ``thick_share``, the share of skeleton
     length wider than 6 px (0.038 +/- 0.012). Branch-level tortuosity
     could not see the problem: the simulated arcades scored on target
     (wide q90 1.082 vs 1.089) while curling back on themselves over 100
     px -- at a 60-px step they turned 30 degrees per step against 15 in
     real eyes. On the spec seed the three terms read 0.58 / 171 px /
     0.119: 129 points of a 159-point score
   - **Finding.** The developmental wave exposes demand only in a thin
     annulus at the front, so a trunk that has caught up with the front
     feels no forward pull, only the sideways pull of the exposed sites
     beside it, and rides along the front: loops concentric with the disc,
     which stay 8-10 px thick because a trunk that circles sheds few
     teeth per unit of radial progress (908 skeleton pixels above 8 px in
     the window against 84 in a real eye)
   - **Model.** Three changes. ``developmental_wave.gate_caliber_max`` 0
     -> 0.006: tips of arcade caliber see the whole hypoxic field and
     head for the far tissue instead of the front's annulus (spec seed
     159 -> 32: alignment 0.58 -> 0.75, reach 171 -> 268 px, thick share
     0.119 -> 0.070, because a trunk that runs sheds teeth and tapers).
     The tenth pass rejected this exemption for breaking the front's
     discipline, under an objective that could not see what it fixes.
     ``boundaries.ArcadeGuidance`` (new component): the astrocyte
     template itself, a force of ``magnitude`` 0.2 pushing tips at or
     above ``min_radius`` 0.006 away from the disc in the plexus plane
     (alignment 0.75 -> 0.83). ``particles.root_radius`` 0.02 -> 0.017:
     trunks start at 8.7 px instead of 10.2 px, against HRF's 7.5-px top
     percentile (thick share 0.087 -> 0.055). Spec seed 159 -> 15
   - Swept and rejected on the spec seed (total vs 159 base): thinner
     roots alone (0.015: 1156, arterial supply 0.34 -- without the gate
     the artery tree starves); a caliber-attenuated hypoxia pull on wide
     tips (``perfusion_demand.caliber_exponent`` 1: 56); more flow per
     comb tooth (``side_branch_flow`` 0.25: 101); guidance alone at 0.3
     (96: straight but thick and short -- the push orients a trunk, the
     demand horizon makes it run); the gate at 0.008 (67, too few tips
     exempt). On top of gate and guidance: ``split_interval`` 21 (44),
     ``radius_taper`` 0.996 (38), ``side_branch_flow`` 0.2 (31),
     ``flow_remodeler.max_radius`` 0.014 (no effect: trunks are never
     adapted). On top of the full candidate: root 0.015 (63, skeleton
     density 2.9%), guidance 0.1 (27) and 0.3 (58, paired perfusion
     0.93: trunks pushed past the tissue)
   - Eight seeds (7/42/909/2024/123456 and 31/77/5150): mean score
     117.3 -> 25.6 (per seed 102/81/108/100/159/95/90/203 ->
     23/19/19/22/15/19/61/27); radial alignment 0.61 -> 0.82, reach 186
     -> 267 px, thick share 0.105 -> 0.063, junction spacing 16.9 -> 20.1
     px (target 20.7); paired perfusion 0.93 -> 0.96 with arterial supply
     >= 0.95 on every seed. Without the guidance force (gate and root
     alone) the mean is the same within noise (27.3) but arterial supply
     falls to 0.88-0.90 on three seeds, which is why the component ships.
     The seed-77 outlier (61) is a thick, far-reaching arcade set (thick
     share 0.095, reach 301 px), not a perfusion failure
   - Held-out seeds 11/202/909/4242 all colonize 100% of the tissue
     (4/4) with arterial supply 0.96/0.98/0.96/1.00 (0.82/0.94/0.99/
     1.00 before). Spec seed: score 159 -> 15.1 under the new objective
     (30.0 -> 15.1 under the old one, which could not see the loops);
     radial alignment 0.58 -> 0.85, reach 171 -> 263 px, thick share
     0.119 -> 0.055, paired perfusion 0.91 -> 0.99 (arterial 0.94 ->
     0.99), skeleton density 3.02% -> 3.69% (target 3.73%), junction
     spacing 16.5 -> 19.8 px (target 20.7), 31k -> 48k frozen particles.
     The residuals are branch tortuosity (1.093 vs 1.079, 3.7 points),
     fractal dimension (1.385 vs 1.349, 2.4) and a wide share below
     target (0.16 vs 0.22, 2.2)

**v0.27.0 - 09/04/26**

 - Paired perfusion: the periphery deficit, and the tip speed behind it

   - **Finding.** Scored as supply *and* drainage, the model perfused
     0.81 of the tissue on the calibration seeds (0.75-0.88) while
     colonizing 0.99 of it. The unsupplied sites were the far
     periphery: 59% of them lay more than 2.5 units from the disc (the
     field ends at 3.2), in the sectors facing away from it, with the
     nearest artery a median 0.27 away -- a tree that stopped short,
     not a hole in the field. Per-type time courses located the cause:
     artery tips died mid-field (median live-tip distance 1.0-1.6 from
     the disc) while vein tips reached the periphery by step 300, and
     both trees ended the run with ~7 live tips, most of them depth-4
     twigs that cannot split again
   - **Model.** Three spec knobs, no new code. ``particles.
     terminal_velocity`` 0.15 -> 0.18: faster tips reach the periphery
     before the mid-field extinction catches them (paired perfusion
     0.81 -> 0.90 on its own), at the price of over-filling the field
     (fractal dimension 1.36 -> 1.39, arcade junction spacing 17.6 ->
     15.0 px against HRF's 20.7). ``path_splitter.split_interval`` 15
     -> 18 and ``side_branch_probability`` 0.65 -> 0.5 pay that back
     with fewer, longer segments and a more open comb (FD 1.366,
     spacing 16.9 px, caliber-profile KS 0.09 -> 0.04) at paired
     perfusion 0.93 (arterial 0.97, venous 0.96). Eight-seed mean score
     103.7 -> 23.3, the paired term 79.8 -> 10.0; the five original
     calibration seeds read 147/126/48/132/67 -> 21/9/13/16/30 and
     three fresh ones 15/10/72 (seed 5150 drains only 0.84 -- the
     deficit can sit on either tree)
   - **Measurement.** ``metrics.arcade_caliber_ratio`` weights each
     depth-0 trunk once inside the zone, as CRAE/CRVE weights each
     vessel, rather than each particle: one sweep seed grew a vein
     trunk that coiled inside the zone while tapering to 1 px, put 147
     of its 526 vein-trunk particles there and read AVR 1.40 -- 210
     points of a 242-point score from one coil
   - Swept and rejected, mean score on 3-5 calibration seeds against a
     base of 103.7: per-type front advance 146 (brakes the healthy
     tree); ``max_depth`` 5 184 (arterial supply 0.66: more twigs
     compete for the same flow); type-scaled steering stiffness 160
     (reverted); ``min_active_tips`` 4 no change; ``cross_type_factor``
     0.1 106; ``influence_radius`` 3.0 83 and ``magnitude`` 0.45 84
     (arteries up, veins down; 80 vs 32 on the spec seed with the
     faster tips); terminal velocity 0.20 41.5 (the FD and spacing
     penalties double); ``side_branch_probability`` 0.4 collapses one
     seed's artery tree (supply 0.67); and a per-type re-sprouting
     watch for the developmental wave (re-sprout the lagging tree on
     its own deficit without holding the front for the other) -- 75.6
     alone, 41.8 with the faster tips against 48.9 without, but it
     collapsed one fresh seed under the shipped comb (154, arterial
     0.75) and seed 909 under all three knobs (104, venous 0.86), so it
     was reverted: which tree lags is a lottery, and forcing sprouts on
     the laggard starves the other
   - Held-out seeds 11/202/909/4242 all colonize 100% of the tissue
     (4/4) with arterial supply 0.82/0.94/0.99/1.00 (0.92/0.99/0.96/
     0.75 before). Spec seed: score 66.6 -> 30.0, paired perfusion
     0.84 -> 0.91 (arterial 0.85 -> 0.94, venous 0.99 -> 0.96),
     caliber-profile KS 0.10 -> 0.04, fractal dimension 1.386 -> 1.355
     (on target), skeleton density 3.33% -> 3.02% (now the largest term
     after perfusion, 5.2 points), branch tortuosity 1.085 -> 1.092
     (3.3 points), pruned particles 12.8k -> 6.2k

**v0.26.0 - 09/04/26**

 - The AVR term, and the starved artery tree behind it

   - **Measurement.** ``artery_vein_caliber_ratio`` is now read on the
     depth-0 arcades within ``metrics.AVR_ZONE`` (0.1-0.5 units) of the
     disc, the analog of the clinical measurement zone, instead of over
     every depth-0 particle. The whole-trunk mean mixed distance bands
     with different type composition and taper, so the ratio read
     0.63-0.85 on trees seeded at 0.67 (a seed lottery of how far each
     tapering trunk ran); in the zone it reads 0.66-0.71 on every seed.
     ``simulation.Geometry`` bundles the containment, perfusion lattice
     and disc position the metrics need (``get_geometry``)
   - **Finding.** The artery tree was 3-6x smaller than the vein tree
     (4,313 vs 26,601 frozen particles on the spec seed) and supplied
     only 59-64% of the tissue on two of three calibration seeds while
     veins drained 99%. Per-lineage bookkeeping traced it to the first
     100 steps: vein trunks spawned 50 new lineages to the arteries' 18,
     mostly second- and third-generation comb teeth. The comb threshold
     (``side_branch_radius`` 0.006) is absolute: vein teeth (0.0106)
     comb for the whole run, artery teeth (0.0071) taper out of comb mode
     after a fifth of it, and the lead then feeds itself (the smaller
     tree's twigs carry the same flow through thinner vessels, so its
     per-tree pruning median is 5x higher and 93% of its terminals fall
     below the pruning threshold against 44% for veins)
   - **Model.** ``path_splitter.type_scaled_comb`` (on): the comb
     threshold is judged on each tree's own caliber scale
     (arteries: ``side_branch_radius x artery_caliber_ratio``), since
     arterioles run at about two thirds of venule caliber at the same
     branching level. Arterial supply 0.64/0.59/0.91 -> 0.85/0.76/0.96 on
     the calibration seeds, artery tree 4.3k -> 10.6k particles, skeleton
     density 3.16 -> 3.33% on the spec seed; held-out seeds 11/202/909/
     4242 all colonize >= 99% of the tissue (4/4) with arterial supply
     0.92/0.99/0.96/0.75. The old-target score is unchanged within noise
     (3-seed mean 17.8 -> 19.7: more 1-px length than the caliber profile
     wants), so the knob is shipped for the arterial supply, not the score
   - **Scoring.** ``perfused_fraction`` now requires both supply and
     drainage: a demand site counts as perfused when it is within
     ``perfusion_radius`` of an artery *and* a vein
     (``metrics.paired_perfused_fraction``). The any-vessel fraction is
     kept as ``colonized_fraction`` (the growth-completeness vital and
     the contact sheet's reliability gate) alongside per-tree coverage.
     This is the honest perfusion and it costs the model on purpose: the
     paired fraction is 0.84 on the spec seed (colonized 0.999), so the
     spec-seed score reads 23.0 -> 66.6, of which 47.8 is this one term
     (the other terms sum to 18.8, down from 23.0 with the AVR term gone)
     -- the largest term in the score again, and it names the next target
   - Swept and rejected on the spec seed, all reverted: a
     supply-gated anastomosis (fuse only where the tip's own tree already
     supplies the tissue; A:V tree ratio 0.16 -> 0.12), a network-wide
     pruning reference instead of per-tree medians (0.16 -> 0.20),
     tissue pressure at the venous end (-0.5: 0.29, -0.7: 0.17),
     equal root calibers (0.14: the feedback, not the caliber, decides
     once a tree is behind), no anastomosis at all (0.64, by letting
     terminal tips run unchecked -- not a fix)

**v0.25.0 - 09/03/26**

 - Thin vessels: the adaptation cap was the 3-px class

   - The re-baselined measurement put 157 of 193 calibration points on
     the caliber profile: half a real eye's fundus-visible skeleton is
     1 px wide, the model's finest superficial vessels were 3 px. The
     cause is ``flow_remodeler.max_adapted_radius``: the shear-adaptation
     set point is each tree's median shear, which the ~20,000 deep-plexus
     capillaries dominate, so nearly every superficial segment reads as
     high-shear and is ground up to the cap -- 52% of superficial length
     sat exactly at 0.006 (3 px). Lowering the cap to 0.003 (1.5 px, a
     precapillary arteriole) moves that mass to the 1-2 px class the real
     eye carries: 3-seed mean 165 -> 27, caliber KS 0.34-0.41 -> 0.03-0.05
     (a real eye scores 0.05), junction spacing back on target, every
     seed still fully perfused
   - ``plexus_layers.dive_radius`` 0.004 -> 0.003: tips of 1.5-2 px are
     precapillary arterioles and venules, which real fundi show in the
     superficial plexus; keeping them there lifts superficial skeleton
     density 2.7% -> 3.1% (HRF 3.7%) and superficial perfusion 76% ->
     95%, 3-seed mean 27 -> 22.5 (23/20/25). Held-out seeds 11/202/909/
     4242 all fully perfused (4/4), superficial skeleton density
     2.3-3.0% -> 3.0-3.4%. Swept and rejected: cap
     0.0035 (32), dive probability 0.02 (37; keeps deep-plexus geometry
     superficial, obtuse share doubles), dive radius 0.0025 (36;
     overshoots the 1-px class), ``max_depth`` 5 (37), side-branch
     probability 0.8 (28)
   - Spec seed: 193.2 -> 23.0. Diameter composition by branch count is
     now 38/44/18% against HRF's 38/40/22%. What remains is spread thin:
     A:V ratio 8.3 (0.81 vs 0.67, a depth-0 arcade-length lottery),
     area density 4.0 (8.8% vs 10.8%), skeleton density 3.4 (3.16% vs
     3.73%), junction spacing 2.2
   - Measurement: ``branch_tortuosity_median`` -> ``branch_tortuosity_mean``
     (clipped at 2). Most branches are a few pixels long, so their
     tortuosity takes a handful of discrete values and the median was a
     lottery among them (13 distinct values across 15 masks; every
     simulation scored the same 1.0706) that no knob could move. The mean
     is continuous (HRF 1.0787 +/- 0.0073; re-derived with
     ``--derive-targets``)

**v0.24.0 - 09/03/26**

 - Re-baseline the V&V measurement: one pipeline for sim and HRF

   - The V&V audit found that much of the calibration score was
     measurement convention rather than model error. Every convention is
     now applied identically to both sources (``vnv.metrics`` module
     docstring): binarization is a majority vote (HRF block-average
     threshold 0.1 -> 0.5; the simulation is drawn by the exact
     pixel-center rule that is the limit of drawing fine and
     majority-downsampling, replacing dilation by a rounded integer
     radius); densities are per imaged pixel (convex hull of the vessel
     pixels) rather than per frame pixel; the skeleton is pruned of spurs
     shorter than a countable branch before anything is counted; branch
     length is the arc length of the pixel chain (diagonal steps sqrt 2,
     so a 45-degree line has tortuosity 1, not 0.71); box counting uses a
     fixed 2-128 px box range instead of one derived from the frame
   - The HRF-derived ``TARGETS`` are re-derived under the new conventions
     and reproducible with ``vnv_calibrate --derive-targets``. The two KS
     targets are now leave-one-out: target and scale are what one real
     eye scores against the other fourteen pooled (log-length 0.052 +/-
     0.028, caliber profile 0.049 +/- 0.029) instead of zero with a
     judgment scale. Branch tortuosity is 1.075 +/- 0.003 (was 1.000, a
     pixel-count artifact); capillary share 0.37 +/- 0.13 (was 0.065 --
     half the HRF skeleton is 1 px wide once thin vessels are no longer
     thickened); wide share 0.22 +/- 0.04 (was 0.345)
   - ``compare.py`` scores through ``calibrate.scoring_stats`` instead of
     its own copy of the statistics; ``image_stats`` measures a real mask
     and a simulated raster with the same function; the perfusion lattice
     parameters are read from the spec (``simulation.get_perfusion_params``)
     rather than hard-coded in three places; the superficial plexus's own
     perfused fraction and the imaged-region share are reported (unscored)
   - Spec-seed score 57.0 -> 193.2, reported plainly and decomposed in the
     roadmap's fifteenth pass: area density 19.1 -> 0.5 and fractal
     dimension 7.9 -> 2.0 were measurement; the skeleton-density gap is
     real (the field of view nets out); and the caliber profile now
     carries 157 points -- the model's finest superficial vessels are 3 px
     wide where half a real eye's skeleton is 1 px, which the old
     threshold had hidden. That is the next modeling target
   - New unit tests for every convention (majority raster widths,
     sub-half-pixel vessels vanishing, majority binarization, hull
     densities, spur pruning, diagonal tortuosity, padding-invariant
     fractal dimension)

**v0.23.0 - 09/03/26**

 - Bug hunt: four PathSplitter/PathFreezer defects fixed

   - A parallel-agent audit of the components, prompted by the seed
     lottery, turned up four bugs of the same species as that one --
     labels that silently change physics, and destructive writes that
     lose vessels without a trace. All four are fixed and unit-tested
     (``tests/test_split_pool.py``), and the healthy model is more
     robust for it: the 3-seed calibration mean holds (67.2 -> 64.8)
     while the per-seed spread tightens from 44-111 to 57-69, every
     seed fully perfused
   - **Free-pool double-claim (critical)**: ``split_paths`` carved the
     free-particle pool once, then let the active-split phase, the
     re-sprout phase, and the two vessel types within the re-sprout loop
     all draw from the same snapshot. When two phases fired in one round
     -- exactly the tree-collapse state the model fights -- they wrote
     new branches to overlapping particle slots, and the later phase's
     write silently overwrote the earlier's, amputating the healthy
     tree's continuations mid-rescue. The pool is now carved
     sequentially across every phase
   - **Shared sprout path id**: ``split_frozen`` never advanced
     ``next_path_id``, so every re-sprout in a round shared one path id
     (and collided with the next split's minor daughter), exempting
     unrelated vessels from ``FrozenRepulsion``'s same-path rule so they
     grew through each other. Each sprout now gets a fresh id
   - **Out-of-plane splits**: the split rotation axis (x_hat x heading)
     tilted out of the plexus plane for tips heading along +/-x,
     throwing daughters into z where the terminal-velocity clamp (z
     weighted heavily by the thin ellipsoid) stalled them; ~43% of tips
     were tilted more than 30 degrees. ``split_axis`` now rotates about
     the plexus normal projected perpendicular to the heading, keeping
     every split in-plane; identical to the old axis for in-plane tips
   - **Dropped OU steering state**: ``PathFreezer`` continuations did
     not carry ``wx/wy/wz``, capping the tortuosity persistence time
     (roadmap idea 7) at the freeze interval regardless of the
     configured value. Continuations now carry it; a no-op at the
     shipped white-noise default
   - Two existing tests held only because of the bugs and are restated:
     ``test_persistence_makes_curvature_coherent`` (was
     ``..._reduces_tortuosity``) asserts what roadmap idea 7 documents --
     persistent steering of fixed spread produces coherent curvature (new
     ``metrics.path_turning_coherence``, the lag-1 autocorrelation of the
     signed turning angle: 0.09 -> 0.54) and *raises* path tortuosity
     (1.010 -> 1.024); the old "straighter" assertion passed only while
     the freezer reset the steering every three steps. The plexus
     stratification test starts its roots on the
     superficial plane, so it measures how well vessels hold their layer
     rather than how fast the spring pulls an off-plane start onto it
     (the in-plane split axis removed the incidental z kicks that had
     been doing that job). In the shipped spec the intermediate and deep
     plexuses stratify better after the fixes (median |z error| 0.022 ->
     0.014 and 0.038 -> 0.034), the superficial slightly worse (0.011 ->
     0.013)

**v0.22.0 - 09/03/26**

 - Root trunks survive their first split: the seed lottery, found and
   fixed

   - Early-trajectory telemetry located the seed lottery: every seed
     has all six root arcades alive at step 25; by step 50 the fit seed
     keeps 5, seed 7 keeps 4, seed 909 keeps 3, and every downstream
     metric scales with the survivors. Two defects, both fixed
   - ``path_splitter.continuation_keeps_path`` (now on):
     ``split_unfrozen`` relabeled the continuing daughter with a new
     ``path_id``, so ``FrozenRepulsion``'s same-path delay exemption no
     longer covered the trail it had just laid -- a dozen springs pushed
     it from behind at full strength and the stacked force crossed the
     extinction threshold at the first split rounds (steps 15-45). The
     continuation now keeps its parent's path. False reproduces legacy
     bit-for-bit, including its quirk that both sisters shared one new id
   - ``PathFreezer.freeze_particles`` froze every tip even when the free
     pool could not supply continuations -- a silent end of the whole
     front. It now tops the pool up and skips the round. This surfaced
     only once trunks survived (faster growth drained the pool): with a
     2000-particle pool seed 7 kept 5 trunks to step 175 and reached
     13.2k particles by step 300 against 3.0k
   - Calibration seeds (incumbent 51.0 / 168.5 / 125.4, mean 115.0):
     44.6 / 111.4 / 45.6, mean 67.2, every seed fully perfused -- the fit
     seed's best score yet, and seed 909, the chronic collapser, now
     matches it. The freezer fix alone (51.0 / 80.5 / 126.7, mean 86.1)
     rescues seed 7 by itself -- the silent mass freeze was hitting it
     under legacy dynamics too -- while the path fix is what lifts the
     fit seed and 909
   - Held-out contact sheet (seeds 11/202/909/4242): perfusion
     0.88 / 1.00 / 0.99 / 0.87 -> 1.00 / 1.00 / 1.00 / 1.00, reliability
     2/4 -> 4/4, mean score 255.3 -> 68.4 (the freezer fix alone reaches
     3/4 at 208.8 / 181.7 / 126.7 / 326.8). For the first time every
     seed the project tracks builds a complete network
   - The bifurcation-angle metric no longer depends on path labels:
     ``true_bifurcations`` recognized a freezer continuation by its
     shared ``path_id``; it now recognizes it by caliber (the one child
     keeping at least 98% of the parent's radius), the structural
     definition, unchanged for legacy trees
   - Falsified on the way: the twelfth pass's radial-emergence
     hypothesis (roots already emerge radially and deterministically),
     and a wider disc circle (radius 0.15 lifts the fit seed to 47.7 but
     leaves seed 909 with no trunks alive by step 100)

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

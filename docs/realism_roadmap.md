# Roadmap: Toward more realistic simulated retinal vasculature

The current model is a *morphological* generator: vessel growth is a persistent
random walk shaped by containment, mutual repulsion, and scheduled splitting.
Real retinal vasculature is shaped by *function* — hypoxia, blood flow, and
metabolic cost. The enhancements below pull the model toward that, roughly in
order of impact per unit of effort. Each maps onto a vivarium component, which
is the point of building this on vivarium in the first place.

After implementing any of them, re-run the V&V harness (see
[Validation & verification](#validation--verification-vv) below) and commit
the overwritten `docs/vnv/` outputs, so the effect of the change is visible
both as an image diff in the pull request and in the quantitative metrics.

## 1. Vessel calibers via Murray's law — IMPLEMENTED

Real branchings satisfy r_parent^k ≈ r_1^k + r_2^k with k ≈ 2.7–3, and branch
*angles* co-vary with radius ratios (Murray's optimality principle: large
daughters deviate little, small side-branches come off near-perpendicular).

*As implemented*: `Particle3D` owns a `radius` column seeded by
`particles.root_radius`; `PathFreezer` continuations inherit the caliber with
a configurable `radius_taper`; `PathSplitter` draws a flow split, assigns
daughter radii by Murray's law (`murray_exponent`, `flow_asymmetry`,
`min_radius`), and derives the branch angles from the radii via the
minimum-work formulas (the configured `split_angle` remains as a fallback for
uncalibered parents); `PathDLA` attachments get a capillary `attach_radius`.
The visualizer, growth GIF, and rasterization all draw true calibers, and the
V&V harness gained two caliber metrics: *vessel area density* and the fitted
*junction exponent* (median k = 3.00 on the generated network; see
`docs/vnv/`). Bifurcation angles now emerge from the radii and concentrate
inside the 60–90° literature band.

## 1b. Caliber-dependent branching cadence — IMPLEMENTED

Real vessel segments keep a roughly constant length-to-diameter ratio: the
arcade trunks run a long way between branch points while capillaries branch
densely. The original splitter checked every tip on the same fixed cadence
(`split_interval` × `split_probability`), so every segment came out about the
same length — visible in the V&V segment-length histogram as an over-tall
spike at the shortest lengths and a deficit through the middle of the
distribution, which was the clearest visual difference from the HRF data.

*As implemented*: `PathSplitter` scales each tip's split probability by
`(min_radius / radius) ** caliber_cadence_exponent` — the configured
`split_probability` applies in full at the capillary caliber floor, and wider
tips split less often in proportion to their caliber (exponent 1 makes the
expected inter-branch distance proportional to diameter, i.e. constant L/D;
0 restores the old caliber-independent cadence). The exponent alone starves
the network (delayed trunk splits compound into far fewer tips), so it is
paired with a faster base cadence: the spec sets `caliber_cadence_exponent`
0.6 with `split_interval` 15 (was 30), so twigs branch more densely than
before while trunks branch less. At the standard 800-step run the
Kolmogorov–Smirnov distance between the sim and HRF log segment-length
distributions fell from 0.099 to 0.053, the median matches HRF exactly
(22 px), the over-tall shortest-length bin dropped from 2.0× to 1.6× the
HRF density, and skeleton density rose from 3.03% to 3.39% (HRF:
3.21 ± 0.29%); as a side benefit the arcade A:V caliber ratio moved from
0.84 to the clinical 0.67. A new tree-based `tree_segment_length` metric in
`metrics.json` tracks inter-branch-point distances in simulation units,
independent of rasterization.

## 2. Growth toward hypoxia (space colonization) — IMPLEMENTED

Real angiogenesis is VEGF chemotaxis: tissue far from any vessel recruits
growth toward it.

*As implemented*: the `PerfusionDemand` force component lays a lattice of
tissue demand sites inside the containment ellipsoid (`site_spacing`); a site
is hypoxic while no frozen vessel lies within `perfusion_radius` (checked
against PathFreezer's KDTree). Following Runions et al.'s space-colonization
rule, each hypoxic site recruits only its *nearest* active growth tip (within
`influence_radius`), and each tip is pulled with strength `magnitude` in the
mean direction of the sites it won — so tips spread apart and fill
unperfused territory instead of clustering, and the force fades away as the
tissue becomes perfused. Keep `magnitude` below `path_extinction`'s
`force_threshold`, which was raised to 1.2 in `model_spec.yaml` so tips
survive long enough to colonize distant tissue. A new `perfused_fraction`
V&V metric measures coverage directly. At the standard 800-step run the
network went from 68.8% to 97.5% tissue perfused, vessel area density from
3.93% to 7.81% (HRF: 11.93 ± 1.03%), skeleton density from 0.96% to 2.48%
(HRF: 3.21%), and skeleton fractal dimension from 1.21 to 1.41 (HRF:
1.35 ± 0.02).

## 3. Paired arterial and venous trees — IMPLEMENTED

The retina has interdigitating artery and vein trees entering at the optic
disc, with an artery:vein caliber ratio around 2:3 (a standard clinical
biomarker).

*As implemented*: a `vessel_type` particle column, seeded as alternating
artery/vein arcades around the disc (artery roots at
`particles.artery_caliber_ratio` × the vein `root_radius`) and inherited down
each tree through freezing, splitting, and DLA attachment. `FrozenRepulsion`
scales cross-tree repulsion by `cross_type_factor` (0.25 in the spec), so
arteries and veins tolerate each other's proximity while avoiding their own
tree — the interdigitating arcade pattern. `PerfusionDemand` became
type-aware: tissue needs both arterial supply and venous drainage, so each
tree is recruited separately and neither can win territory for both; and
`PathSplitter` re-sprouts any tree whose active tips all died (angiogenic
sprouting), closing the ratchet that let one tree take over. The visualizer
and growth GIF color the trees like a fundus photo, and the V&V harness
gained per-tree coverage metrics and the arcade A:V caliber ratio; the
tree-based metrics now measure only *true* bifurcations (same-path
continuation children are excluded).

## 4. Anastomosis: close the loops — IMPLEMENTED

Trees don't perfuse; circuits do. Real retinal capillaries form a mesh
connecting the arterial and venous trees. `PathDLA` already gestures at this
(stray particles attaching to frozen vessels). Make it explicit: when an
active tip of one tree comes within a capillary radius of the *other* tree's
terminal segments, freeze a connecting segment. This produces a perfusable
graph and enables flow-based remodeling (below).

*As implemented*: the `PathAnastomosis` component fuses an active
capillary-caliber tip (radius ≤ `max_tip_radius`) onto the other tree when
it comes within `capture_radius` of an opposite-type frozen segment of
capillary caliber (≤ `max_target_radius`) — capillaries join capillaries,
tips don't fuse into trunks. The tip freezes with `probability` per step
once in range, and the join is recorded in a new `anastomosis_id` particle
column, so the result is a real graph edge: `tree_edges`, the rasterizer,
the visualizer, and the growth GIF all draw the bridges (violet in the
renders), and idea 5's flow solve has a closed circuit to work on. A new
`graph_cycles` V&V metric (E − N + C over parent and anastomosis edges)
counts the independent loops, which is exactly the number the tree-only
model could never move off zero. At the standard 800-step run the spec
settings (`capture_radius` 0.03, caliber caps 0.004, `probability` 0.5)
produce 40 anastomoses forming 37 independent loops, with the idea-1b
headline metrics essentially unchanged (skeleton density 3.33%, 98.7%
perfused, KS to the HRF length distribution 0.068) — wider capture radii
or caliber caps were swept and actually produce *fewer* loops, because
early fusions kill tips and shrink the network before it interdigitates.

## 5. Flow-based remodeling and pruning — IMPLEMENTED

The strongest realism lever in the literature (Pries & Secomb's adaptation
model): over-produce segments, solve Poiseuille flow on the frozen graph (a
sparse linear solve — cheap at this scale), then prune segments below a
wall-shear threshold and thicken high-flow ones. Real development does exactly
this. The schema is already waiting: the `unfreeze_time` column exists and is
unused — pruning is what it was born for. A `FlowRemodeler` component running
on a slow cadence (every N steps, like `PathSplitter`) fits the event system
naturally.

*As implemented*: the `FlowRemodeler` component (see
[docs/poiseuille_flow.md](poiseuille_flow.md) for the physics from first
principles) treats every frozen segment and anastomosis bridge as a resistor
of conductance r⁴/L, fixes artery roots at `artery_pressure` and vein roots
at `vein_pressure`, gives every node a small `leak_fraction` conductance to
tissue (so dead-end twigs carry a trickle instead of exactly zero), and
solves Kirchhoff's current law as one sparse linear system per
`remodel_interval`. Each segment's wall shear |Q|/r³ then drives two
adaptations: *pruning* — terminal segments (degree-1 graph ends, so the
network is never cut mid-branch) with shear below
`shear_threshold_fraction` × their own tree's median are recycled back into
the free particle pool, stamping the long-dormant `unfreeze_time` column —
and *caliber adaptation* — every segment's radius drifts by
`adaptation_rate` toward the caliber that would put its shear at its own
tree's median, so high-flow trunks thicken and low-flow twigs thin. Shear
targets are per tree because arteries genuinely run at higher wall shear
than veins; a global target thickens arteries until the clinical A:V
caliber ratio inverts. The anastomosis loops from
idea 4 carry real artery-to-vein flow and are naturally protected by their
high shear; the low-shear dead-end clutter that the diameter-stratified V&V
flagged (~20× too many capillary-scale branches) is exactly what gets eaten
back, tip by tip, pass after pass.

At the standard 800-step run the spec settings (`shear_threshold_fraction`
0.5, `adaptation_rate` 0.15) prune 3,281 segments over the run while the
network keeps growing (7,137 frozen segments remain, 98.9% perfused): the
capillary-scale share of skeleton branches falls from 66% to 27% (HRF:
6.5%), vessel area density jumps from 7.2% to 15.1% as vessels thicken
(HRF: 11.9 ± 1.0% — we overshoot now instead of undershooting), skeleton
fractal dimension moves from 1.49 to 1.43 (HRF: 1.35), and the arcade A:V
caliber ratio holds at 0.66 (clinical ~0.67). Costs, stated plainly: the
aggregate segment-length median drifts from 22 to 29 px (HRF: 22) as the
shortest clutter branches are removed, and the junction exponent is no
longer the imposed k = 3.00 — adapted radii give an *emergent* median
k ≈ 2.2 (n = 41), which is in fact closer to what fundus measurements
report for human retina (means around 2.2–2.6) than the theoretical 3.
Two design points mattered: shear targets are per tree (see above), and
milder pruning or pruning without adaptation were swept and do worse —
the two mechanisms work as a pair.

## 6. Layered plexuses — IMPLEMENTED

Retinal vasculature is stratified — superficial, intermediate, and deep
capillary plexuses connected by short diving vessels. The thin ellipsoid is
currently one layer. Quantize z into 2–3 sub-layers with type-dependent layer
preference and occasional vertical connections, and the model can generate
synthetic OCT-angiography stacks, not just fundus-style projections.

*As implemented*: a `layer_id` particle column names each vessel's home
plexus, inherited through freezing, splitting, and DLA like `vessel_type`.
The `PlexusLayers` component holds every active tip near its layer's
z-plane (`layer_z`, superficial first) with a *damped* Hookean spring —
undamped, tips ping-pong through the plane, saturate the terminal-velocity
clamp, stall, and go extinct — and each step a capillary-caliber tip
(radius ≤ `dive_radius`) dives one layer deeper with `dive_probability`;
the frozen trail left during the transit is the diving vessel. Arteries,
veins, and wide arterioles never dive, so the deeper plexuses are
capillary-only, matching the anatomy. Flattening a 3D slab onto planes
concentrates the frozen-repulsion crowding and stacks the layer force onto
the extinction budget, so `path_extinction.force_threshold` rose to 2.0 and
the anastomosis `capture_radius` widened to 0.045 (wide enough to reach
across a plexus gap). At the standard 800-step run the three plexuses hold
5,079 / 931 / 1,624 segments connected by 173 diving vessels, with the
established headline metrics intact: 97.9% perfused, skeleton density
2.69%, vessel area density 14.2%, fractal dimension 1.43, arcade A:V ratio
0.68, and 41 anastomoses forming 39 loops (up from 21). A new
`docs/vnv/plexus.png` renders OCTA-style en-face slabs per layer plus an
x–z cross-section, and `metrics.json` gains a `plexus_layers` block
(per-layer counts, calibers, z-adherence, diving-vessel count).

## 7. Smooth, controllable tortuosity — IMPLEMENTED

The white-noise velocity jitter produces wiggle at the step scale; real
vessels wander smoothly, and tortuosity *changes* are diagnostic (retinopathy
of prematurity, diabetic retinopathy). Replace the random velocity term with
an Ornstein–Uhlenbeck (autocorrelated) process: one interpretable knob — a
persistence time — that directly controls a clinically meaningful biomarker.

*As implemented*: `Particle3D` carries an OU steering state (`wx/wy/wz`
columns) updated as an AR(1) with correlation time
`particles.noise_persistence_time` (in days) and the same stationary
spread as the old uniform kick (sd = `overall_max_velocity_change`/√3), so
the process degenerates exactly to white noise as the persistence time
approaches one step — and `noise_persistence_time: 0` keeps the legacy
white-noise path bit-for-bit.

What the 800-step sweep taught us, stated plainly: the healthy model's
tortuosity was *already* at the target (branch median 1.004 vs. HRF 1.000)
— white-noise kicks average out at the branch scale, so the healthy spec
keeps `noise_persistence_time: 0` and its calibration unchanged. Turning
persistence on is the *pathology* dial: correlated steering produces
coherent curvature, raising path tortuosity (1.011 → 1.049 at 0.25 days)
— and because persistent tips commit into repulsion walls instead of
diffusing around them, the phenotype arrives as a package with capillary
dropout and hypoperfusion (97.9% → 75–90% perfused), reminiscent of
diabetic retinopathy's tortuous, dropout-ridden vasculature. Amplitude
compensation (holding the integrated wander constant, sd ∝ 1/√τ) does not
decouple the two: the network's force ecology is genuinely tuned around
diffusive tips. A clean single-variable tortuosity dial for *healthy*
networks would need the retuning pass of idea 8; as a disease knob it
works today.

## 8. Calibrate to real statistics — then break them deliberately — HEALTHY FIT IMPLEMENTED

Pick validation targets from public datasets (HRF, DRIVE, STARE fundus vessel
segmentations; OCTA-500 for depth-resolved data) — fractal dimension (~1.7 for
healthy retina), vessel density, branch-angle and segment-length
distributions, tortuosity, foveal avascular zone area, junction exponents —
and fit the configuration parameters against them. Once the healthy model is
calibrated, disease phenotypes become parameter perturbations: capillary
dropout and neovascular tufts (diabetic retinopathy), a peripheral avascular
front (ROP), collateral formation after simulated vein occlusion. If the goal
is synthetic training data, the *pathological* variety is where the value is.

*As implemented (healthy fit)*: the `vnv_calibrate` CLI formalizes the
objective the per-feature sweeps had been eyeballing — eight validation
targets (skeleton density, area density, fractal dimension, branch
tortuosity from the HRF across-mask means and sds; KS distance to the HRF
segment-length distribution; capillary branch share; clinical A:V ratio
0.67; full perfusion, one-sided) each scored as a squared z-like deviation,
summed into a single calibration score where 0 is perfect and each unit is
one squared sd-equivalent of miss. `metrics.json` now carries the score of
every standard run (`calibration` block, plus the total in the figure
headline) so it is tracked across model versions like every other metric.
The fit itself is coordinate descent over the eight knobs the sweeps showed
the metrics respond to, one 800-step run per evaluation, every evaluation
logged.

The first fit (22 evaluations) moved three knobs —
`caliber_cadence_exponent` 0.6 → 0.75, `adaptation_rate` 0.15 → 0.10,
`perfusion_demand.magnitude` 0.3 → 0.35 — and cut the calibration score
from 38.7 to 24.7 (−36%): KS distance to the HRF length distribution
0.171 → 0.074, capillary branch share 22.5% → 19.1%, fractal dimension
1.43 → 1.41 (HRF 1.35), area density 14.2% → 13.6% (HRF 11.9%), with the
A:V ratio (0.66) and tortuosity (1.003) holding; skeleton density gave
back 2.69% → 2.43% (HRF 3.21%) as the search traded it against the
distributional terms. The search surface is genuinely coupled and sharp:
`split_interval` 18 collapses the network outright (score ≈ 4,500), and
the old cadence exponent 0.6 scores 340 under the *new* incumbent — the
knobs cannot be tuned one feature at a time, which is exactly why the
harness exists. Robustness: the fitted knobs beat the previous spec on
*every* seed tested (38.7 → 24.7 on the fit seed; 269 → 157 and 46 → 45
on two held-out seeds), but the seed-to-seed spread is larger than the
within-seed improvement — whether a run catches a good growth trajectory
dominates — so a natural next step for the harness is a multi-seed
objective (average the score over 2–3 seeds per evaluation, at
proportional compute cost). Disease phenotypes also remain future work:
the ingredients (the OU tortuosity dial, per-type perfusion, pruning)
are in place.

*Second fit (diameter composition)*: visual review after the first fit
showed the simulated calibers splitting into thin-and-thick with too few
mid-caliber vessels, so the composition became an explicit target.
Binning skeleton branches by their EDT-recovered diameter (≤2 px,
2–4 px, >4 px at the standard raster scale), HRF is 6/59/35%
(capillary/mid/wide) while the first fit produced 19/15/66% — the mid
bin, which dominates real networks, was the *smallest* simulated
stratum. The comparison figure now carries a diameter histogram and
grouped composition bars (by branch count and by skeleton length), and
the wide share joined the calibration targets. A composition-aware
coordinate-descent pass found no single-knob improvement; the per-target
logs exposed why — shear adaptation *polarizes* calibers away from the
median (below-median twigs thin toward the floor while above-median
vessels thicken), actively hollowing out the mid stratum — so the
remodeler gained an `adaptation_deadband` (segments within a factor of
their tree's median shear don't adapt; the default 1.0 preserves the
previous behavior). The winning move in the follow-up sweeps was gentler
remodeling rather than the deadband itself: `adaptation_rate`
0.10 → 0.05, `shear_threshold_fraction` 0.5 → 0.65, anastomosis
`capture_radius` 0.045 → 0.035, moving the composition to 21/21/58% and
the score 40.2 → 34.8 (these totals include the new wide-share target,
so they are not comparable to the 24.7 above). The three knobs interact:
on the fit seed, reverting any one of them is worse than the trio
(40.4 / 85.8 / 46.0). Robustness is mixed in the same way as the first
fit: the trio also wins on a held-out healthy seed (54.8 → 40.5, with
the composition moving the same direction), but on the known-degenerate
seed 42 — where both configs produce a broken sparse network — the
stiffer pruning threshold and slower adaptation amplify the degeneracy
(165 → 326, mostly unperfused tissue). The fix for that is the
multi-seed objective above, not a different point fit. The residual gap is
structural rather than a tuning miss: Murray-law bifurcations plus the
geometric freeze taper transit the mid-caliber band in a couple of
branch generations, so real mid-bin mass needs either longer
mid-caliber runs between branchings or a caliber-dependent taper —
mechanism work for a future pass, ideally calibrated against the
per-plexus ROSE data once the Zenodo request is granted.

*Third pass (arcade geometry)*: at the zoomed-out scale the wide vessels
meandered and curled, unlike the smooth HRF arcades. Two mechanisms were
responsible, and both were fixed at the mechanism level rather than by
tuning. First, steering was caliber-blind — an arcade tip took the same
random kicks as a capillary tip, so its heading random-walked — fixed by
`noise_caliber_reference/_exponent`: tips wider than the capillary
reference have their random steering attenuated by
`(reference/radius)^exponent`, so arcades hold their heading while
capillaries wander (real large vessels are stiff; sprouting is a
capillary-tip behavior). Second, and less obvious: the flow remodeler
had an *anastomosis-shortcut runaway* — an artery→vein capillary bridge
sees an enormous pressure drop through a tiny radius, its shear
explodes, adaptation thickens it, conductance grows as r⁴, it draws yet
more flow — a positive feedback that promoted curly capillary paths all
the way to arcade caliber (on the fit seed, only 376 of 2,516
wide-caliber particles were true arcades, and 191 sat in the
capillary-only deep plexuses). The new `max_adapted_radius` saturates
shear-driven *thickening* at venule caliber while segments born wider
can still taper. Alongside the mechanisms, the HRF comparison itself was
corrected to be fundus-faithful: fundus photographs do not see the deep
capillary-only plexuses (OCTA does — see `plexus.png`), so the
comparison now rasterizes the superficial layer only, and with the
promotion runaway gone the missing wide mass turned out to be *real
arcade mass* — seeding 6 root trunks instead of 4 closes most of the
superficial density gap and finally makes the network panel read like a
fundus image: smooth arcades sweeping from the disc into a fine mesh.
The honest cost: the calibration score prefers the old meandering
config (39.0 vs 55.0 on the corrected measurement) because skeleton-
branch tortuosity is chopped at junctions and cannot see the
long-wavelength curl — the fat promoted vessels bought area density and
wide-share cheaply. The `wide_tortuosity_q90` target (HRF 1.11,
one-sided) guards the worst of it, but a proper vessel-tracking
curliness metric (merge branches through junctions along the widest
continuation, then measure direction drift per unit arc) is the missing
instrument; until it exists, the steering exponent and adaptation cap
stay pinned during refits instead of being left to the score.

*Fourth pass (comb-like side branching)*: even with smooth arcades, the
branches off them were far too sparse — real arcades are *monopodial*: a
trunk that keeps nearly its own caliber and sheds small side branches at
short, comb-like intervals, where the splitter only did near-symmetric
dichotomous forks and the caliber cadence made wide tips branch rarely.
The pattern now has its own instrument: `wide_junction_spacing_px`, the
skeleton distance between branch points along wide (>4 px) vessels
(junction clusters counted once, measured identically on sim rasters and
HRF masks). HRF carries a branch point every 22.7 px of wide skeleton
(sd 1.8); the previous model managed one every 38.6 px — an 80-point
score term that quantifies exactly what the eye saw. The splitter's new
comb mode (`side_branch_flow` / `side_branch_radius` /
`side_branch_probability`, off by default) makes parents above the
reference caliber emit side branches at their own cadence and strongly
asymmetrically: the trunk keeps ~96% of its caliber and the tooth takes
the Murray caliber for a ~10% flow fraction, leaving near-perpendicular
on a random side — both angles fall out of the minimum-work bifurcation
relations already in the code. The instructive failure: the first comb
collapsed the network outright, because teeth spawn inside the frozen-
repulsion field of their own trunk (interaction radius 0.15 units
= 38 px — wider than the comb spacing itself) and the stacked forces
drove every tip extinct; the repulsion radius, it turns out, was itself
an anti-realism constraint, and dropping it to 0.12 both rescues the
comb and lets vessels pack at real densities. The fit (emission 0.65,
spacing 26.0 px vs HRF 22.7) scores 47.3 against the previous spec's
135.7 on the corrected objective, and the emission rate was chosen for
robustness, not just the fit seed: 0.8 packs tighter (19.9 px, 50.5)
but tips one held-out seed into the crowding-extinction cascade, while
0.65 holds on both. Superficial skeleton and area density remain the
weakest terms — the fundus-visible network is still sparser than HRF —
and the operating point is sensitive: small pushes toward more
superficial mass (lower dive probability, weaker repulsion) destabilize
the growth front.

*Fifth pass (multi-seed objective)*: the harness investment every
earlier pass asked for. `vnv_calibrate --seeds 123456,7,42` makes the
objective the *mean* score across the listed seeds (one simulation per
seed per evaluation, per-seed breakdowns logged), so a config that
collapses on any seed loses to one that is merely mediocre everywhere.
Re-judging the comb-era candidates under it was clarifying: every
candidate has exactly one weak seed (the growth front either catches a
good trajectory or it doesn't), so single-seed fits had been silently
trading robustness for fit-seed polish. The best 3-seed mean moved one
knob — `dive_probability` 0.05 → 0.04, keeping slightly more of the
capillary bed superficial — improving the mean 114.1 → 70.6 (per-seed
47/74/221 → 42/121/49) and the fit seed itself. Candidates that looked
tempting on one seed (a stiffer extinction threshold, gentler tooth
emission, a wider capture radius, and their combinations) all lost on
the mean. The remaining seed-to-seed spread (42 vs 121) is honest
variance in whether the arcades colonize evenly; averaging over more
seeds per evaluation (at proportional compute) or seed-averaged
coordinate descent with a larger budget is the obvious continuation
when compute allows.

*Sixth pass (bifurcation angles)*: the tree-based angle histogram was
bimodal — a peak near 75° plus a second mode at 100–150° that healthy
retinas don't show. Decomposing every junction by provenance located the
second mode precisely: it is almost entirely *deep-plexus capillary
junctions* (layer-2 obtuse share 0.80 vs 0.08 superficial; capillary
parents median 116° vs arcade/comb parents 76° with zero obtuse share) —
sibling tips pulled toward opposite hypoxic voids splay into the
T-shaped junctions that real deep capillary meshes genuinely have, but
that fundus photographs cannot see. The fundus-visible geometry was
already unimodal on the Murray optimum. So the fix is the same
measurement principle as the superficial raster: the figure's angle
panel and Murray-exponent inset now measure the superficial tree only
(the all-layer summary stays in `metrics.json`), and two literature
targets pin it — `bifurcation_angle_median` 77±5° and one-sided
`bifurcation_obtuse_share` (>100°, target 0.05). The current spec sits
at median 75–79° and obtuse 0.09–0.14 across the three calibration
seeds. A mechanism dial was also built and honestly rejected:
`perfusion_demand.caliber_exponent` attenuates hypoxia chemotaxis on
wide tips (biologically plausible — VEGF acts on capillary sprouts),
but the multi-seed objective vetoed it decisively (3-seed mean
72.3 → 2563 at exponent 1.0): the pull on wide tips is what drives the
arcades outward, and removing it stalls colonization. The knob remains
(default off) as a candidate disease dial.

*Seventh pass (the length-weighted caliber profile)*: the diameter
strata (≤2 / 2–4 / >4 px) were always a coarse lens, so the harness
gained the binning-free version: `skeleton_pixel_diameters` measures the
local (2×EDT) diameter at every skeleton pixel, making its distribution
exactly *skeleton length by width*; the figure's diameter panel now
plots this caliber profile and `ks_caliber_profile` (KS against the
pooled HRF profile) joined the targets. The profile localized the
residual mismatch to two features the strata could not see: a pile-up
at ~2.8 px — vessels parked exactly at `max_adapted_radius` — and a 3×
deficit in the 5.5–7 px band that real fundi carry. The tooth autopsy
explained both: vein comb teeth are born at 4.8 px *wanting to thicken*
(their shear sits above the tree median), but the cap only lets them
shrink, so 76% of depth-1 particles had been ground to ≤3.2 px with the
median exactly at the cap. The multi-seed fit that fixed it:
`adaptation_deadband` 1.0 → 2.0, so moderate-shear segments keep their
born caliber instead of being dragged to the cap, plus
`side_branch_flow` 0.1 → 0.15, so vein teeth are born at ~5.4 px inside
the missing band (3-seed mean 79.7 → 67.4, caliber KS better on every
seed, obtuse-angle share halved). Raising the cap itself re-triggers
the shortcut runaway (means 244/161, with or without the deadband) and
stays rejected. One deliberate trade is on record: deadband-alone edges
the aggregate mean (64.2) by winning back seed-42 density terms, but
loses the caliber profile and angle terms this pass targets — the
thicker-teeth config ships because matching length×width was the goal.

*Eighth pass (sub-teeth, and the percolation problem)*: visual review
still read the mid vessels as too long and meandering between branch
points — comb teeth ran ~90 px between splits (the cadence damping)
where HRF mid vessels branch every 30–50 px, and artery teeth sat just
below `side_branch_radius`, so they never combed.
`side_branch_radius` 0.008 → 0.006 extends the comb one caliber class
down (teeth grow sub-teeth), and on the fit seed this produces the most
fundus-like texture yet, with near-target skeleton density (2.9%) and
full perfusion. The trade is on record deliberately: the 3-seed mean
regresses (67.4 → 165) because one seed in three stalls mid-growth
under the denser branching. Nine stabilizer configurations were swept
under the multi-seed objective — higher and intermediate extinction
thresholds, tighter repulsion, deep-only anastomosis, a self-healing
growth front (`min_active_tips`), developmental gating of the comb
(`side_branch_start_time`), and a doubled particle pool (falsified
bit-for-bit: idle wanderers exert no forces and draws are keyed, so the
run is unchanged) — and each rescued one seed while sinking another.
The conclusion worth recording: the growth front is percolation-like
(early crowding near the disc decides whether a seed's network fills),
so branching density and robustness are coupled through a critical
point, and no per-knob tuning decouples them. Candidate real fixes for
a future pass: growth in waves (freeze the front's advance rate),
crowding-aware tooth emission (skip a tooth when local frozen density
is high), or accepting per-seed retries as the cost of a
near-critical growth model. The three stabilizer knobs stay in the
code, default-legacy and unit-tested, for that work and for disease
phenotypes.

*Ninth pass (growth reliability — a negative result, caught by the new
held-out artifact)*: a field study on four fresh seeds (11/202/909/4242,
none used in calibration) confirmed the percolation prediction — two of
four networks stalled mid-growth. Three candidate fixes from the eighth
pass's list were implemented and unit-tested: a **crowding gate**
(`max_crowding` within `crowding_radius`) that skips dichotomous splits
into saturated space — the immediate-death daughters whose extinction
cascade collapses marginal seeds — with comb teeth exempt, since a
trunk's own frozen trail dominates its neighbor count and the ungated
version wiped out the combs (arcade spacing 55 px vs the 22.7 px
target); an **established self-healing front** (`min_active_tips` gated
by `resprout_established_size`) that re-sprouts a thinning front from
frozen vessels once a tree has a real front to lose (an ungated floor
tops up at the start, crowds the disc, and cost the fit seed 3x); and
**balanced arterial inflow** (`flow_remodeler.balanced_arterial_inflow`:
inject equal flow at each artery root instead of fixing root pressures,
so no arcade starves its siblings — the rich-get-richer failure). An
ten-config sweep found a pair (gate 30 + floor 8 established at 200)
that improved the 3-seed calibration mean 581.4 → 431.8 with no
calibration seed worse — and the new held-out contact sheet then showed
reliability falling from 2/4 seeds ≥95% perfused to 0/4. Attribution
probes pinned each harm: the tip floor alone drops seed 11 from 0.996
to 0.76 perfused and seed 4242 from 0.928 to 0.54 (re-sprouted tips
near frozen mass seed the very cascades they were meant to heal); the
gate at 30 drops seed 202 from 0.987 to 0.79, and it cannot be tuned
around — healthy interiors run 30–44 frozen neighbors within 0.06
(p90–max), pathological pile-ups no higher, so any gate low enough to
fire also caps normal filling, and a gate of 60 never fires at all
(bit-for-bit identical runs). Balance was rejected separately: on
collapse-prone seeds it feeds flow into stalled fragments (gate+floor
mean 504.9 → 673.6 with it on). Everything therefore ships default-off,
and the conclusion sharpens the eighth pass's: per-seed outcomes near
the critical point are chaotic, and per-tip local rules reshuffle which
seeds fill rather than making every seed viable. The still-promising
directions are the non-local ones — growth in waves (control the
front's advance rate globally) or accepting per-seed retries. What the
pass does ship enabled is measurement: `vnv_contact_sheet` renders the
four held-out seeds beside HRF masks with a [STALLED] stamp and records
the reliability fraction in `contact_sheet.json`, so this problem is
tracked across versions instead of rediscovered — it is exactly what
caught the non-generalization above.

*Tenth pass (the developmental wave — growth as a closed loop)*: the
ninth pass ended pointing at non-local mechanisms, and the first one
tried — growth in waves, the way real retinal vascularization spreads
from the optic disc behind a hypoxia front — is the first mechanism in
the project's history to improve every calibration seed at once. The
`DevelopmentalWave` component owns a front radius around the disc:
`PerfusionDemand` exposes only demand sites within `radius + lookahead`,
so every tip chases the same expanding ring (and unserved holes behind
the front stay visible until filled); the front advances only while the
tissue behind it is served, so a stalling seed concentrates demand at
the stall instead of failing silently; and a front held too long
re-sprouts the stalled tree from the frozen vessels nearest its
unserved sites — the targeting the ninth pass's global tip floor
lacked. The control variable is the goal itself (perfusion behind the
front), not a per-tip local proxy, and that is exactly why it escapes
the ninth pass's no-go result. One discovery mattered most in tuning:
with per-type advance the front brakes on the artery tree's ~0.85
mid-field service ceiling (arteries alone genuinely cannot cover
tissue at the perfusion radius), and a healthy seed held at the brake
is ground down by pruning — `advance_rule: combined` (any-vessel
service) removes the brake and restores the fit seed bit-for-near
(47.0 → 51.0 at full perfusion). The swept landscape is a textbook
single-knob tension: threshold 0.75 outruns weak seeds (909: 376), 0.95
makes strong seeds wait (seed 7: 432), 0.85 balances (3-seed mean
581.4 → 115.0, with chronic collapser 909 at 125.4 and 0.99 perfusion —
a seed no per-tip mechanism ever moved above 0.45). A wide-tip
exemption (arcades seeing the whole field) tested worse: it breaks the
front discipline. On the held-out contact sheet the catastrophic
collapse mode is eliminated — min perfusion 0.40 → 0.87, mean
0.83 → 0.94, mean score 574.7 → 255.3 — with one honest trade: strict
≥95% reliability stays 2/4, because two previously-fast seeds settle at
a stable ~0.87–0.88-perfused equilibrium instead of finishing. A
long-horizon probe (seed 11 to step 1200) shows the plateau is real,
not pacing: from step 800 on, shear pruning consumes what the held
front builds, and the system balances just below the bar. So the wave
converts legacy's bimodal outcome (0.99 or 0.40, a coin flip per seed)
into a compressed 0.87–1.00 band — the floor rises dramatically, two
ceilings come down slightly. The natural follow-up is to close the
remaining gap on the pruning side: FlowRemodeler currently prunes
low-shear terminals everywhere, including the young sprouts at a held
front that the wave is trying to grow; sparing the front region from
pruning (or coupling `shear_threshold_fraction` to front distance)
should let the equilibrium finish. The wave ships enabled.

*Eleventh pass (the prune grace period — a negative result that
corrects the tenth pass's diagnosis)*: the follow-up was built as the
simplest possible lever, `flow_remodeler.prune_grace_days` — terminal
segments frozen less than this long ago are not pruned, so a young
sprout gets time to connect and earn its shear. Age-based sparing was
chosen over a position band around the front because it needs no
cross-component coupling and also protects the wave's behind-front
re-sprouts (one knob, four lines, unit-tested). The sweep on the
calibration seeds rejected every value: 2.5 days is worse on all three
(fit 51.0 → 60.1, seed 7 168.5 → 241.7, seed 909 125.4 → 154.8); 5 and
10 days rescue seed 7 spectacularly (→ 80.3 / 86.0 at full perfusion,
22–26k particles) but halve the fit seed's network (→ 121.5 / 126.6,
perfusion 0.87–0.89, capillary share 0.10 → 0.35) and collapse seed 909
(125.4 → 580.0 / 641.4, perfusion 0.59 / 0.55); the 3-seed means are
152.2 / 260.6 / 284.7 against the incumbent 115.0. The mechanism
corrects the tenth pass's reading: pruning was not only consuming the
held front's sprouts, it was **decluttering** — a low-shear dead-end
stub that is spared stays frozen and acts as a `FrozenRepulsion` source
that chokes the tips around it, so whether spared stubs mature into
vessels (seed 7) or into clutter (the fit seed) is seed-specific. That
is the ninth pass's per-seed reshuffle signature again, now on a
remodeling knob rather than a growth knob, and the shipping state stays
the wave's 0.87–1.00 band. The two 0.88 seeds therefore sit at a
tip-extinction equilibrium at the held front, not a pruning one; the
next lever to try is on the extinction side (a front-local relief of
repulsion or of the force threshold for sprouts born from a
`resprout_toward_stall` call), and it should be gated by the held-out
contact sheet like everything since the ninth pass. The pass also
consolidated the stacked increments without changing behavior
(bit-for-bit over 150 steps): one eligibility filter and one commit path
in `PathSplitter`, one nearest-vessel query shared by `PerfusionDemand`
and the wave, one disc-distance helper, and one seed/build helper pair
for the V&V CLIs.

*Twelfth pass (hypoxic survival relief — the extinction-side lever,
another negative result, and what three of them together say)*: the
eleventh pass's corrected diagnosis — a tip-extinction equilibrium at
the held front — was attacked at the extinction side with the most
biologically grounded lever available: hypoxia-driven VEGF is a
survival signal for endothelial tip cells as well as a chemoattractant,
so tips sitting in unserved tissue should tolerate more force before
regressing. `PathExtinction` now reads its per-tip threshold from a
`particle.extinction_threshold` value pipeline and `PerfusionDemand`
raises it by `survival_factor` for tips within `perfusion_radius` of a
hypoxic site of their type (the pipeline wiring avoids an import cycle
and keeps the mechanism with the demand field it depends on; 1.0 is
bit-for-bit legacy). The mechanism it targets is real: a sprout from
`resprout_toward_stall` gets a new `path_id`, so `FrozenRepulsion`'s
own-path `delay` exemption does not cover its frontier parent, and the
frontier's summed repulsion pushes it over the threshold before it can
move. The sweep rejected it anyway: the fit seed is bit-identical across
factors 1.5/2/3 (51.0 → 61.8, still fully perfused — the relief is
effectively binary, the same tips are saved at every factor), seed 7
improves enormously (168.5 → 70.1 / 67.0 / 67.0 at full perfusion), and
seed 909 collapses at every factor (125.4 → 967.2 / 884.2 / 884.2,
perfusion 0.99 → 0.44 / 0.47, network 11098 → 3240–4035 particles);
3-seed means 366.3 / 337.7 / 337.7 against the incumbent 115.0.

That is now three independent levers on three subsystems — growth
(ninth pass: crowding gate, tip floor), remodeling (eleventh: prune
grace), extinction (twelfth: survival relief) — each of which rescues
seed 7 and collapses seed 909, or the reverse, and only the wave (a
non-local controller on the goal itself) lifted every seed. The
consistent reading: the per-seed outcome is decided early, by the
trajectory the six root arcades take out of the disc, and post-hoc
per-tip relief anywhere downstream only reshuffles which trajectory
wins. The initial root headings are the one seed-dependent input that
has never been examined: `initial_velocity_range` draws each root's
direction at random, so some seeds start with arcades aimed inward or
bunched, while real arcades emerge from the disc radially in four
quadrants. A radially-outward emergence with a small angular jitter is
the next candidate — cheap, developmentally faithful, and the first
lever that acts on the cause rather than the symptom — and it should
be gated on the held-out contact sheet like everything since the ninth
pass.

*Thirteenth pass (the seed lottery, found and fixed)*: the twelfth
pass's hypothesis died on contact with the code — `initialize_circle_positions`
already sends every root out radially, evenly spaced and deterministically;
`initial_velocity_range` only randomizes the idle pool, which exerts no
forces. So the divergence was measured instead of guessed: an
early-trajectory telemetry (active tips, frozen count, trunks alive,
extinctions every 25 steps) showed all three calibration seeds identical
at step 25 with six root arcades alive, and by step 50 the fit seed
keeping 5, seed 7 keeping 4, seed 909 keeping 3 — with every downstream
number scaling with the survivors. The seed lottery is trunk survival
through the first split rounds (steps 15–45). A wider disc circle (0.15,
putting roots outside each other's repulsion radius) was falsified next:
the fit seed improves to 47.7 but seed 909 has no trunks alive by step
100. The actual killer is in `split_unfrozen`: both daughters of a split
received a *new* `path_id`, including the continuing trunk, so
`FrozenRepulsion`'s same-path delay exemption no longer covered the
trail the trunk had just laid — a dozen Hookean springs within 0.12
pushed it from behind at full strength, and the stacked force crossed
the extinction threshold. Which trunks split first, and so which died,
was the comb-emission draw: the whole lottery. The fix is one line of
semantics (`continuation_keeps_path`: the continuing daughter keeps its
parent's path; legacy behavior, sisters sharing one new id and all, is
preserved bit-for-bit when off) — and the moment trunks survived, a
second defect surfaced: `PathFreezer.freeze_particles` froze every tip
even when the free pool could not supply continuations, silently ending
the whole front (seed 7: 33 active tips → 2 in one freeze round, no
extinction marks). Faster growth drained the pool first, which is why
this never showed under legacy dynamics and why the ninth pass's
pool-size probe was bit-identical. It now tops the pool up and skips the
round. Together (incumbent 51.0 / 168.5 / 125.4, mean 115.0): **44.6 /
111.4 / 45.6, mean 67.2**, every seed fully perfused — the fit seed's
best score yet, and seed 909, unmovable above 0.45 perfusion by any
per-tip lever, now matches it with a 41k-particle network. The freezer
fix alone (51.0 / 80.5 / 126.7, mean 86.1) rescues seed 7 by itself —
the silent mass freeze was hitting it under legacy dynamics too — while
the path fix is what lifts the fit seed and 909. On the held-out contact
sheet (seeds 11/202/909/4242) perfusion goes 0.88 / 1.00 / 0.99 / 0.87 →
1.00 / 1.00 / 1.00 / 1.00, reliability 2/4 → 4/4, mean score 255.3 →
68.4 (the freezer fix alone reaches 3/4, at 208.8 / 181.7 / 126.7 /
326.8): for the first time every seed the project tracks builds a complete
network, and the percolation problem the eighth pass named as central is
closed. One measurement had to follow the semantics:
`true_bifurcations` recognized freezer continuations by shared
`path_id` and went NaN under the fix; it now recognizes them by caliber
(the one child keeping ≥98% of the parent's radius — split daughters are
born at Murray calibers below that), the structural definition,
unchanged for legacy trees. The lesson the three negative passes were
pointing at was right — the outcome is decided early — but the cause
was not a parameter at all; it was a labeling convention that turned a
vessel's own trail into its repulsor.

*Fourteenth pass (a bug hunt, and what it says about the measurements)*:
finding the seed lottery in a labeling convention rather than a
parameter prompted a systematic audit — parallel agents over the growth
components, the force/flow/plexus components, and the V&V measurement
code, each hunting the same species (labels that silently change
physics, destructive writes, stale reads, sim-vs-real measurement
mismatch). The growth audit found four more bugs, all fixed and
unit-tested: a **free-pool double-claim** where the active-split and
re-sprout phases wrote new branches to overlapping particle slots and
the later silently overwrote the earlier (the same destructive-write
species as the freezer bug, firing in the tree-collapse state); a
**shared sprout path id** (re-sprouts never advanced the id counter, so
unrelated vessels were mutually exempt from repulsion — the same
label-decides-physics species as the split relabeling); **out-of-plane
splits** (the rotation axis tilted daughters into z for x-heading tips,
where the terminal-velocity clamp stalled them, ~43% of tips affected);
and a **dropped OU state** at freezer handoff (capping the tortuosity
persistence time). The healthy model came out more robust — the 3-seed
mean held (67.2 → 64.8) while the per-seed spread tightened (44–111 →
57–69), held-out reliability stayed 4/4. The force/flow/plexus audit
cleared the suspected prune-orphans-a-tip hazard (zero occurrences over
15,871 prunes) and left minor latent notes (a dormant force-cache
staleness, a soft-core disc repulsion whose magnitude rises with
distance). The measurement audit is the consequential one: it found
that **roughly half the calibration score is measurement artifact, not
model error** — the simulation's density is computed over the
ellipsoid's bounding box (whose corners are ~22% empty by construction)
while the HRF frame is nearly filled, and the HRF masks are binarized
with a threshold that thickens real vessels ~28% while the simulation is
rasterized crisp, so the caliber and composition terms agree partly by
cancellation. Measuring both sides on a common field of view and a
common binarization would collapse the area-density gap from 38% to 8%
and re-baseline most of the distributional targets. That re-baselining
is a deliberate, judgment-heavy change to the calibration foundation and
is tracked as its own follow-up rather than folded into the bug fixes.

*Fifteenth pass (the measurement re-baseline)*: every convention the
audit flagged is now applied identically to both sources (the V&V
section below lists them), and the HRF-derived targets were re-derived
under the new conventions with `vnv_calibrate --derive-targets`. The
score of the spec seed moves 57.0 → 193.2, and the decomposition says
where the audit was right and where it was wrong. Right: the
**area-density term was measurement** (19.1 → 0.5; the 0.1 binarization
threshold had thickened HRF vessels so the target was 11.9% instead of
10.8%, and the simulation's dilation-drawn calibers were rounded to odd
widths), the **fractal term flipped sign** (sim above target by 3 sd →
below by 1.5 sd once both sides count the same box range — a 256-px box
that only the simulation's square frame admitted was inflating its
dimension), and the tortuosity target of 1.000 was an artifact of
counting pixels rather than arc length (real branches sit at 1.075, and
the model at 1.071 turns out slightly *too straight*). Wrong: the field
of view nets out (HRF masks fill 83% of their frame, the ellipse 80%), so
the skeleton-density gap is real — the model lays 29% less centerline
per imaged pixel than a healthy eye (2.66% vs 3.73 ± 0.31%). And the
re-baseline revealed the term that now dominates everything: the
**length-weighted caliber profile** (KS 0.41 against a real eye's 0.05 ±
0.03, 156 of the 193 points). Half of a real eye's fundus-visible
skeleton runs in vessels one pixel wide at this scale (~15–30 µm —
precapillary arterioles and venules, which the HRF labelers drew and
which the old threshold had inflated to 3 px, hiding the mismatch);
the model's finest superficial vessels are its side branches at radius
0.006 (3 px), and only 5% of its branches reach the 1-px class against
38% in HRF. The composition and density gaps are the same fact seen
three ways. This is the modeling target for the next pass, and it is a
caliber-generation question — thinner side branches, a deeper terminal
generation at the caliber floor, or a shear-adaptation floor below the
current `min_radius` — not a growth-reliability one; every seed already
perfuses. Two more things the re-baseline surfaced, unscored: the
superficial plexus alone perfuses 76% of the tissue (the deep plexuses
carry the rest, consistent with anatomy), and the deep plexus sits a
median 0.034 from its plane against a 0.04 layer gap, i.e. it barely
stratifies — a plexus-spring question for a later pass.

*Sixteenth pass (thin vessels: the adaptation cap was the 3-px class)*:
the autopsy of the caliber profile took one histogram. Binning the
superficial frozen particles by radius, 52% sat *exactly* at
`flow_remodeler.max_adapted_radius` = 0.006 (3 px) — not near it, at it.
The mechanism: shear adaptation drifts every segment toward its tree's
*median* shear, and the median is set by the ~20,000 deep-plexus
capillaries, most of which the same adaptation has already thinned to
the 0.001 floor (68% of deep vessels are sub-pixel). Against that
median nearly every superficial segment reads as high-shear, so the
deadband never spares it and it is ground up to the cap; the cap, not
Murray's law, was setting the caliber of the fundus-visible fine
vessels. The seventh pass had seen this pile-up and treated it as a
feature to fill the 5–7 px band — which the old 0.1 binarization
threshold had manufactured by thickening HRF's thin vessels. Lowering
the cap to 0.003 (1.5 px, a precapillary arteriole) moves that mass into
the 1–2 px class: 3-seed mean 165 → 27, caliber KS 0.34–0.41 → 0.03–0.05
(a real eye scores 0.05 ± 0.03 against the others), junction spacing
14 → 19 px (target 21), every seed fully perfused. The second knob
follows from the same anatomy: `plexus_layers.dive_radius` 0.004 →
0.003, so tips of 1.5–2 px — precapillary arterioles and venules, which
fundus photographs show in the superficial plexus — stay superficial
instead of diving; superficial skeleton density 2.7% → 3.1% (HRF 3.7%),
superficial perfusion 76% → 95%, mean 27 → 22.5 (23/20/25, the tightest
spread the project has had). Swept and rejected on the calibration
seeds: cap 0.0035 (32), dive probability 0.02 (37 — it keeps deep-plexus
mesh geometry superficial and the obtuse-angle share doubles), dive
radius 0.0025 (36 — overshoots the 1-px class, KS back to 0.13),
`max_depth` 5 (37), side-branch probability 0.8 (28). On the spec seed
the score is 193 → 23 and the diameter composition by branch count is
38/44/18% against HRF's 38/40/22%. The remainder is spread thin and
partly noise: the A:V ratio (0.81 vs 0.67) is a depth-0 lottery of how
far each tapering trunk runs; area density 8.8% vs 10.8% and skeleton
density 3.16% vs 3.73% say the superficial network is still ~15% short
of length. One measurement fix rode along: the branch-tortuosity
*median* was quantized — most branches are a few pixels long, their
arc/chord takes a handful of discrete values, and every simulation
scored the identical 1.0706 (13 distinct values across the 15 masks) —
so the target is now the clipped mean (HRF 1.0787 ± 0.0073), which is
continuous and which knobs can move. The lesson is the same one the
measurement pass taught from the other side: a caliber class that
appears in *both* the model and the reference at the same pixel width
is not evidence the model is right, when a cap on one side and a
threshold on the other are what put it there.

*Seventeenth pass (the AVR term, and the starved artery tree behind
it)*: the A:V term (8.3 points, 0.81 against a configured 0.67) was a
measurement artifact of a familiar kind — the estimator averaged every
depth-0 particle, so it mixed distance bands with different type
composition and taper, and read anywhere from 0.63 to 0.85 depending on
how far each trunk happened to run. Read in a zone near the disc, as
the clinical CRAE/CRVE is, it is 0.66–0.71 on every seed. Behind it,
though, the per-type diagnostics exposed a real asymmetry: the artery
tree was 3–6× smaller than the vein tree and supplied only 59–64% of
the tissue on two of three calibration seeds while veins drained 99%.
Six single-seed counterfactuals and a per-lineage bookkeeping run
located it. Not the caliber thresholds (equal root calibers: still
7× asymmetric, because a tree that falls behind stays behind), not the
pruning reference (a network-wide median instead of per-tree: 0.16 →
0.20), not the tissue pressure, not anastomosis consuming the fusing
tips (a gate that fused only in tissue the tip's own tree already
supplied made it worse). The bookkeeping found the asymmetry in the
first 100 steps — vein trunks spawned 50 lineages to the arteries' 18,
mostly second- and third-generation comb teeth — and the comb rule is
why: `side_branch_radius` is one absolute number, vein teeth (0.0106)
sit above it for the whole run while artery teeth (0.0071) taper below
it after a fifth of the run, so the vein tree gets a comb of combs and
the artery tree one generation less. From there the lead feeds itself:
both trees carry the same flow, the smaller one through thinner and
fewer vessels, so its per-tree pruning median runs 5× higher and 93% of
its terminal twigs fall below the threshold against 44% for veins.
`type_scaled_comb` judges the comb threshold on each tree's own caliber
scale (arteries at `side_branch_radius × artery_caliber_ratio`, the
observation that arterioles run at about two thirds of venule caliber
at every level): arterial supply 0.64/0.59/0.91 → 0.85/0.76/0.96 on the
calibration seeds, the artery tree 4.3k → 10.6k particles, superficial
skeleton density 3.16% → 3.33%. The score itself moves little under
the old targets (mean 17.8 → 19.7: more 1-px length than the caliber
profile wants), which exposed the last blind spot in the objective —
`perfused_fraction` counted any vessel, so a tissue fed by veins alone
scored as perfused. It now requires both supply and drainage (a site
within reach of an artery *and* a vein), the any-vessel fraction stays
as the colonization vital and the contact sheet's reliability gate, and
the paired fraction (0.84 on the spec seed, against 0.999 colonized) is
deliberately the largest term in the score — 47.8 of 66.6, the other
terms summing to 18.8 where the whole score was 23.0 before: the network
colonizes its territory but does not yet perfuse a sixth of it, and
that, not density, is the next target.

## Validation & verification (V&V)

Idea 8 is where every other idea gets measured, so the repo carries a V&V
harness under `vivarium_eye_vessels.vnv`:

- `vnv_growth_gif` renders an animated GIF of the vessel formation process
  from any model specification, for qualitative before/after comparison.
- `vnv_compare` runs a model specification to completion, computes network
  metrics (box-counting fractal dimension, skeleton vessel density, segment
  length distribution — aggregate and stratified by local vessel diameter,
  which is recovered from the distance transform along the skeleton so real
  masks without explicit calibers are measured the same way as the sim —
  per-branch tortuosity, bifurcation angles), computes the same image-based
  metrics on expert-labeled vessel masks from the public HRF dataset
  (downloaded on first use), and writes a side-by-side diagnostic figure
  plus a `metrics.json` for quantitative tracking across model versions.
  The simulation panel shows a fundus-sized central window of the raster at
  the HRF working image's pixel extent (calibers are matched in pixels, so
  this is the like-for-like magnification), with the full field as an inset.
  `metrics.json` also records wall-clock runtime (setup, simulation,
  steps/second) so speed regressions surface alongside the network metrics;
  compare runtimes only within the same machine, and expect them to scale
  with network size.
- `vnv_contact_sheet` runs the model on held-out seeds (defaults: 11, 202,
  909, 4242 — none used in calibration) and renders each superficial network
  beside expert-labeled HRF masks in the same style, stamping every
  simulation panel with its reliability vitals (perfused fraction, skeleton
  density) and a `[STALLED]` flag when perfusion misses 95%. One healthy
  retina looks much like another, so a realistic model must produce a usable
  network on *every* seed, not just the calibration seeds;
  `contact_sheet.json` records the per-seed vitals and the seed-reliability
  fraction so this is tracked across versions like every other metric.

The current model's outputs live in the single `docs/vnv/` folder
(`growth.gif`, `comparison.png`, `metrics.json`, `contact_sheet.png`,
`contact_sheet.json`). After implementing a change, regenerate them in place
with:

```bash
vnv_growth_gif src/vivarium_eye_vessels/model_specifications/model_spec.yaml
vnv_compare src/vivarium_eye_vessels/model_specifications/model_spec.yaml
vnv_contact_sheet src/vivarium_eye_vessels/model_specifications/model_spec.yaml
```

(both default to writing into `docs/vnv/`; keep `--steps` at its default of
800 so runs are comparable). Committing the overwritten files makes each
model change's before/after visible as a side-by-side image diff in the pull
request, and prior versions remain retrievable from git history.

*Measurement conventions* (the module docstring of `vnv.metrics` is the
reference; every convention is applied identically to the simulation and to
the HRF masks, which is what makes the comparison apples to apples):

- Binarization is a majority vote — a pixel is vessel when more than half of
  it is covered. HRF masks are block-averaged to the common 1024-px raster
  and thresholded at 0.5; simulated segments are drawn by the exact
  pixel-center rule (a center within the segment's radius) that is the limit
  of drawing at infinite resolution and majority-downsampling. Vessels
  narrower than half a pixel vanish from both, as they do from a fundus
  photograph, and vessels near one pixel wide fragment in both.
- Densities are per *imaged* pixel — the convex hull of the vessel pixels —
  not per frame pixel, since a fundus camera's circular field and the
  model's elliptical territory both leave empty corners (both fill ~83% of
  their frame, as it happens).
- The skeleton is pruned of spurs shorter than a countable branch (5 px)
  before anything is counted, so ragged edges do not manufacture branch
  points; branch length is the arc length of the pixel chain (diagonal steps
  count √2), so a straight 45° line has tortuosity 1. Branch tortuosity is
  scored by its mean (clipped at 2), not its median: most branches are a
  few pixels long, so the median lands on one of a handful of discrete
  values.
- Box counting uses a fixed range of box sizes (2–128 px) rather than one
  derived from the frame size.
- The HRF-derived targets are the across-mask mean and sd of each metric
  under these conventions; the two KS targets are leave-one-out — what one
  real eye scores against the other fourteen pooled. `vnv_calibrate
  --derive-targets` reprints them, so a convention change re-derives the
  targets in one step.
- The A:V caliber ratio is read on the depth-0 arcades within a fixed zone
  of the disc (`metrics.AVR_ZONE`), as the clinical CRAE/CRVE is, not over
  each trunk's whole tapering run. Perfusion requires both supply and
  drainage: a demand site is perfused when an artery *and* a vein lie
  within `perfusion_radius`; the any-vessel fraction is reported as the
  colonized fraction and gates the contact sheet's seed reliability.

*Comparability caveats that remain*: real fundus masks are 2D projections of
a curved surface while the sim is projected from a thin ellipsoid; bifurcation
angles are measured in 3D on the tree (the plexus is nearly planar, so they
agree with the fundus projection the literature reports); and the sim's
spatial scale is arbitrary, pinned to HRF only through the pixel scale of the
calibers. The harness is designed to show the *direction and size of
improvement*, not to claim the current model matches reality.

## References

- Murray CD. The physiological principle of minimum work. PNAS 1926.
- Pries AR, Secomb TW. Structural adaptation of microvascular networks. Am J
  Physiol 1998.
- Runions A et al. Modeling and visualization of leaf venation patterns
  (space colonization). SIGGRAPH 2005.
- Budai A et al. Robust vessel segmentation in fundus images (HRF dataset).
  Int J Biomed Imaging 2013.
- Masters BR. Fractal analysis of the vascular tree in the human retina. Annu
  Rev Biomed Eng 2004.

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
  `metrics.json` also records wall-clock runtime (setup, simulation,
  steps/second) so speed regressions surface alongside the network metrics;
  compare runtimes only within the same machine, and expect them to scale
  with network size.

The current model's outputs live in the single `docs/vnv/` folder
(`growth.gif`, `comparison.png`, `metrics.json`). After implementing a
change, regenerate them in place with:

```bash
vnv_growth_gif src/vivarium_eye_vessels/model_specifications/model_spec.yaml
vnv_compare src/vivarium_eye_vessels/model_specifications/model_spec.yaml
```

(both default to writing into `docs/vnv/`; keep `--steps` at its default of
800 so runs are comparable). Committing the overwritten files makes each
model change's before/after visible as a side-by-side image diff in the pull
request, and prior versions remain retrievable from git history.

*Comparability caveats*: topological metrics are computed on skeletons for
both sim and real masks so they don't depend on calibers; real fundus masks
are 2D projections of a curved surface while the sim is projected from a thin
ellipsoid; and the sim's spatial scale is arbitrary, so scale-dependent
metrics (density, segment lengths) are normalized to the field of view. The
harness is designed to show the *direction and size of improvement*, not to
claim the current model matches reality — the remaining density gap versus
HRF is capillary-bed structure, which idea 4 (anastomosis) and finer-scale
branching are what close.

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

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

## 4. Anastomosis: close the loops

Trees don't perfuse; circuits do. Real retinal capillaries form a mesh
connecting the arterial and venous trees. `PathDLA` already gestures at this
(stray particles attaching to frozen vessels). Make it explicit: when an
active tip of one tree comes within a capillary radius of the *other* tree's
terminal segments, freeze a connecting segment. This produces a perfusable
graph and enables flow-based remodeling (below).

## 5. Flow-based remodeling and pruning

The strongest realism lever in the literature (Pries & Secomb's adaptation
model): over-produce segments, solve Poiseuille flow on the frozen graph (a
sparse linear solve — cheap at this scale), then prune segments below a
wall-shear threshold and thicken high-flow ones. Real development does exactly
this. The schema is already waiting: the `unfreeze_time` column exists and is
unused — pruning is what it was born for. A `FlowRemodeler` component running
on a slow cadence (every N steps, like `PathSplitter`) fits the event system
naturally.

## 6. Layered plexuses

Retinal vasculature is stratified — superficial, intermediate, and deep
capillary plexuses connected by short diving vessels. The thin ellipsoid is
currently one layer. Quantize z into 2–3 sub-layers with type-dependent layer
preference and occasional vertical connections, and the model can generate
synthetic OCT-angiography stacks, not just fundus-style projections.

## 7. Smooth, controllable tortuosity

The white-noise velocity jitter produces wiggle at the step scale; real
vessels wander smoothly, and tortuosity *changes* are diagnostic (retinopathy
of prematurity, diabetic retinopathy). Replace the random velocity term with
an Ornstein–Uhlenbeck (autocorrelated) process: one interpretable knob — a
persistence time — that directly controls a clinically meaningful biomarker.

## 8. Calibrate to real statistics — then break them deliberately

Pick validation targets from public datasets (HRF, DRIVE, STARE fundus vessel
segmentations; OCTA-500 for depth-resolved data) — fractal dimension (~1.7 for
healthy retina), vessel density, branch-angle and segment-length
distributions, tortuosity, foveal avascular zone area, junction exponents —
and fit the configuration parameters against them. Once the healthy model is
calibrated, disease phenotypes become parameter perturbations: capillary
dropout and neovascular tufts (diabetic retinopathy), a peripheral avascular
front (ROP), collateral formation after simulated vein occlusion. If the goal
is synthetic training data, the *pathological* variety is where the value is.

## Validation & verification (V&V)

Idea 8 is where every other idea gets measured, so the repo carries a V&V
harness under `vivarium_eye_vessels.vnv`:

- `vnv_growth_gif` renders an animated GIF of the vessel formation process
  from any model specification, for qualitative before/after comparison.
- `vnv_compare` runs a model specification to completion, computes network
  metrics (box-counting fractal dimension, skeleton vessel density, segment
  length distribution, per-branch tortuosity, bifurcation angles), computes
  the same image-based metrics on expert-labeled vessel masks from the public
  HRF dataset (downloaded on first use), and writes a side-by-side diagnostic
  figure plus a `metrics.json` for quantitative tracking across model
  versions.

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

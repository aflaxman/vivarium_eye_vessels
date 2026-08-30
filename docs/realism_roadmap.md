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

## 2. Growth toward hypoxia (space colonization)

Real angiogenesis is VEGF chemotaxis: tissue far from any vessel recruits
growth toward it. The elegant version here is a "space colonization" force —
the model already keeps a KDTree of frozen particles, so a `PerfusionDemand`
component can attract each active tip toward the direction of greatest
distance from existing vessels (sampled on a coarse lattice inside the
ellipsoid). Tips then fan out naturally, fill territory evenly, and terminate
when no unperfused tissue remains — a principled, emergent replacement for
`PathExtinction`'s force threshold and much of the repulsion tuning.

## 3. Paired arterial and venous trees

The retina has interdigitating artery and vein trees entering at the optic
disc, with an artery:vein caliber ratio around 2:3 (a standard clinical
biomarker). Add a `vessel_type` to paths: same-type repulsion strong,
cross-type repulsion weak. The characteristic alternating A/V arcade pattern
falls out, and renders can color the trees like a fundus photo.

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
claim the current model matches reality — expect the current model to remain
measurably far from HRF in fractal dimension and density (ideas 2 and 4 above
are what close that gap).

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

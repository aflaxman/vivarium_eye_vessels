"""Calibrate the healthy model against real-data targets (roadmap idea 8).

The objective scores a simulated network against the validation targets in
``TARGETS``, drawn from the HRF dataset (means and spreads of the
image-based metrics) and from clinical literature (the arcade A:V caliber
ratio, full perfusion):
each component is a squared z-like deviation, and the total is their sum,
so a score of 0 means every metric sits on its target and each unit is one
squared standard-deviation-equivalent of miss.

Usage::

    vnv_calibrate src/vivarium_eye_vessels/model_specifications/model_spec.yaml \
        --budget 24 --log calibration_log.json

The HRF-derived entries of ``TARGETS`` are reproducible: ``vnv_calibrate
--derive-targets`` measures the 15 masks with the current
:mod:`~vivarium_eye_vessels.vnv.metrics` conventions and prints the
across-mask mean and sd of every image target, plus the leave-one-out
spread of the two KS statistics (each eye against the pooled others),
which is what a real eye scores and therefore the target and scale for
the simulation. Re-run it whenever a measurement convention changes.

The search is coordinate descent over ``SEARCH_SPACE``: one knob at a time,
try the candidate values around the current setting, keep the best, move to
the next knob, and repeat until the budget runs out or a full pass makes no
improvement. Every evaluation (config, per-metric scores, total) is
appended to the log as it happens, so partial runs are still informative.
Runs are deterministic given the spec's seed; use ``--seed`` to check a
candidate's robustness on other seeds, or ``--seeds 123456,7,42`` to make
the objective the mean score across several seeds (one simulation per seed
per evaluation), so the fit cannot win by overfitting a lucky growth
trajectory.
"""

import copy
import itertools
import json
import time
from pathlib import Path

import click
import numpy as np
import yaml
from scipy.stats import ks_2samp

from vivarium_eye_vessels.vnv import metrics, reference_data, simulation

# (target, scale): score component is ((value - target) / scale) ** 2.
# HRF-derived targets use the across-mask mean and sd (``--derive-targets``
# reprints them under the current measurement conventions); clinical targets
# use literature values with a judgment scale; "one_sided" components only
# score deviations in the bad direction.
TARGETS = {
    "skeleton_density": {"target": 0.0373, "scale": 0.0031},
    "area_density": {"target": 0.1080, "scale": 0.0098},
    "fractal_dimension": {"target": 1.3489, "scale": 0.0229},
    # Arc length over chord of skeleton branches, mean (clipped at 2). Real
    # branches bend a little between junctions; a chord-straight network is
    # as wrong as a meandering one. The mean, not the median: most branches
    # are a few pixels long, so their tortuosity takes a handful of discrete
    # values and the median is a lottery among them (13 distinct values
    # across 15 masks) that no model knob can move
    "branch_tortuosity_mean": {"target": 1.0787, "scale": 0.0073},
    # KS statistics: the target and scale are what a real eye scores against
    # the pooled other eyes (leave-one-out mean and sd), one-sided because
    # matching the pool more closely than a real eye does costs nothing
    "ks_log_length": {"target": 0.0522, "scale": 0.0277, "one_sided": "above"},
    # Diameter composition by branch count. Half the HRF skeleton runs in
    # 1-px vessels, so the capillary share is large but varies widely
    # between eyes (a knife-edge on whether a thin vessel is 1 or 2 px);
    # shares sum to 1, so scoring the capillary and wide shares also pins
    # the mid (2-4 px) share
    "capillary_share": {"target": 0.3698, "scale": 0.1311},
    "wide_share": {"target": 0.2223, "scale": 0.0388},
    # Real arcades run straight: the 90th-percentile tortuosity of wide
    # (>4 px) branches. One-sided: only meandering wide vessels are penalized
    "wide_tortuosity_q90": {"target": 1.0892, "scale": 0.0070, "one_sided": "above"},
    # Comb-like side branching: real arcades carry a branch point every
    # ~21 px of wide (>4 px) skeleton (connected junction clusters counted
    # once, spurs pruned)
    "wide_junction_spacing_px": {"target": 20.73, "scale": 1.97},
    # Length-weighted caliber profile: KS between the per-skeleton-pixel
    # diameter distributions (sim superficial raster vs pooled HRF) — the
    # binning-free version of the composition targets, matching
    # "length x width" across the whole caliber range
    "ks_caliber_profile": {"target": 0.0491, "scale": 0.0286, "one_sided": "above"},
    # Fundus-visible (superficial-tree) bifurcation geometry: healthy
    # arteriolar branch angles are unimodal around ~75-84 degrees; obtuse
    # (>100 degree) junctions are rare. Judgment scales from the literature
    "bifurcation_angle_median": {"target": 77.0, "scale": 5.0},
    "bifurcation_obtuse_share": {"target": 0.05, "scale": 0.05, "one_sided": "above"},
    # Clinical AVR, read on the depth-0 arcades within metrics.AVR_ZONE of the
    # disc (the measurement zone), not over each trunk's whole tapering run
    "artery_vein_caliber_ratio": {"target": 0.67, "scale": 0.05},
    # Tissue within perfusion_radius of both an artery and a vein: a bed
    # with supply but no drainage (or the reverse) is colonized, not perfused
    "perfused_fraction": {"target": 0.98, "scale": 0.02, "one_sided": "below"},
}

# Candidate values per knob, current spec setting included. Chosen from the
# per-feature sweeps: these are the knobs the headline metrics respond to.
SEARCH_SPACE = {
    ("particles", "noise_caliber_exponent"): [0.0, 0.75, 1.5],
    ("path_splitter", "split_interval"): [12, 15, 18],
    ("path_splitter", "caliber_cadence_exponent"): [0.45, 0.6, 0.75],
    ("path_splitter", "side_branch_flow"): [0.06, 0.1, 0.15],
    ("path_splitter", "side_branch_probability"): [0.5, 0.65, 0.8],
    ("path_splitter", "side_branch_radius"): [0.005, 0.006, 0.008],
    ("frozen_repulsion", "interaction_radius"): [0.08, 0.1, 0.15],
    ("path_freezer", "radius_taper"): [0.994, 0.996, 0.998],
    ("flow_remodeler", "shear_threshold_fraction"): [0.35, 0.5, 0.65],
    ("flow_remodeler", "adaptation_rate"): [0.05, 0.10, 0.15],
    ("flow_remodeler", "adaptation_deadband"): [1.0, 2.0, 4.0],
    ("flow_remodeler", "max_radius"): [0.012, 0.016, 0.02],
    ("flow_remodeler", "max_adapted_radius"): [0.006, 0.008, 0.010],
    ("plexus_layers", "dive_probability"): [0.035, 0.04, 0.05],
    ("path_anastomosis", "capture_radius"): [0.035, 0.045, 0.055],
    ("frozen_repulsion", "spring_constant"): [1.25, 1.5, 1.75],
    ("perfusion_demand", "magnitude"): [0.25, 0.3, 0.35],
    ("perfusion_demand", "caliber_exponent"): [0.0, 0.5, 1.0],
}


def calibration_score(stats: dict) -> dict:
    """Per-target squared deviations and their total, NaN-tolerant.

    Missing or NaN stats score as a 5-sigma miss so broken runs never win.
    """
    scores = {}
    for name, spec in TARGETS.items():
        value = stats.get(name)
        if value is None or not np.isfinite(value):
            scores[name] = 25.0
            continue
        deviation = (value - spec["target"]) / spec["scale"]
        one_sided = spec.get("one_sided")
        if one_sided == "above":
            deviation = max(deviation, 0.0)
        elif one_sided == "below":
            deviation = min(deviation, 0.0)
        scores[name] = float(deviation**2)
    scores["total"] = float(sum(scores.values()))
    return scores


HRF_IMAGE_TARGETS = (
    "skeleton_density",
    "area_density",
    "fractal_dimension",
    "branch_tortuosity_mean",
    "capillary_share",
    "wide_share",
    "wide_tortuosity_q90",
    "wide_junction_spacing_px",
)


def image_stats(image: dict, references: dict | None = None) -> dict:
    """The scored image statistics of one binary vessel image.

    ``image`` is a :func:`~vivarium_eye_vessels.vnv.metrics.image_metrics`
    result; with ``references`` (pooled HRF ``lengths`` and
    ``pixel_diameters``) the two KS statistics are included. The same
    function measures a simulated raster and a real mask, which is what
    makes the HRF targets and the simulation's statistics comparable.
    """
    from vivarium_eye_vessels.vnv.compare import stratify_by_diameter

    lengths = np.asarray(image["branch_length_px"], dtype=float)
    tortuosity = np.asarray(image["branch_tortuosity"], dtype=float)
    diameter = np.asarray(image["branch_diameter_px"], dtype=float)
    pixel_diameters = np.asarray(image["pixel_diameter_px"], dtype=float)
    strata = stratify_by_diameter(lengths, diameter)
    wide_tortuosity = tortuosity[diameter > 4.0]
    stats = {
        "skeleton_density": image["skeleton_density"],
        "area_density": image["area_density"],
        "fractal_dimension": image["fractal_dimension"],
        "branch_tortuosity_mean": (
            float(np.clip(tortuosity, 1.0, 2.0).mean()) if len(tortuosity) else float("nan")
        ),
        "capillary_share": len(strata["diameter_le_2px"]) / max(len(lengths), 1),
        "wide_share": len(strata["diameter_gt_4px"]) / max(len(lengths), 1),
        "wide_tortuosity_q90": (
            float(np.quantile(wide_tortuosity, 0.9)) if len(wide_tortuosity) else float("nan")
        ),
        "wide_junction_spacing_px": image["wide_junction_spacing_px"],
    }
    if references is not None:
        stats["ks_log_length"] = (
            float(ks_2samp(np.log10(lengths), np.log10(references["lengths"])).statistic)
            if len(lengths) and len(references["lengths"])
            else float("nan")
        )
        stats["ks_caliber_profile"] = (
            float(ks_2samp(pixel_diameters, references["pixel_diameters"]).statistic)
            if len(pixel_diameters) and len(references["pixel_diameters"])
            else float("nan")
        )
    return stats


def hrf_references() -> dict:
    """HRF reference data: pooled distributions plus per-mask image metrics.

    ``lengths`` and ``pixel_diameters`` pool the 15 masks for the KS
    targets; ``per_mask`` keeps each mask's full
    :func:`~vivarium_eye_vessels.vnv.metrics.image_metrics` result (with
    its ``file`` name) for the per-mask statistics and figures.
    """
    per_mask = []
    for path in reference_data.fetch_hrf_masks():
        binary = metrics.binarize_mask(reference_data.load_mask(path))
        image = metrics.image_metrics(binary)
        image["file"] = path.name
        per_mask.append(image)
    return {
        "lengths": np.concatenate(
            [np.asarray(m["branch_length_px"], dtype=float) for m in per_mask]
        ),
        "pixel_diameters": np.concatenate(
            [np.asarray(m["pixel_diameter_px"], dtype=float) for m in per_mask]
        ),
        "per_mask": per_mask,
    }


def derive_hrf_targets(references: dict) -> dict:
    """Across-mask mean and sd of every HRF-derived target, under current conventions.

    The KS entries are leave-one-out: each mask's statistic against the
    other masks pooled, so the target is what a real eye scores.
    """
    per_mask = references["per_mask"]
    table = {
        name: np.array([image_stats(m)[name] for m in per_mask]) for name in HRF_IMAGE_TARGETS
    }
    loo_lengths, loo_diameters = [], []
    for i, mask in enumerate(per_mask):
        others = {
            key: np.concatenate(
                [np.asarray(m[field], dtype=float) for j, m in enumerate(per_mask) if j != i]
            )
            for key, field in (
                ("lengths", "branch_length_px"),
                ("pixel_diameters", "pixel_diameter_px"),
            )
        }
        stats = image_stats(mask, others)
        loo_lengths.append(stats["ks_log_length"])
        loo_diameters.append(stats["ks_caliber_profile"])
    table["ks_log_length"] = np.array(loo_lengths)
    table["ks_caliber_profile"] = np.array(loo_diameters)
    return {
        name: {"target": float(np.mean(values)), "scale": float(np.std(values))}
        for name, values in table.items()
    }


def scoring_stats(pop, edges, geometry: simulation.Geometry, references: dict) -> dict:
    """The scored summary statistics for one finished simulation.

    ``geometry`` carries the model's containment, perfusion lattice, and
    disc position (:func:`~vivarium_eye_vessels.vnv.simulation.get_geometry`).
    Besides the scored entries, the result carries unscored vitals: frozen
    count, loops, the imaged-region share of the frame, the any-vessel
    colonized fraction (whole network and superficial plexus alone), and
    the coverage of each tree.
    """
    # Fundus photographs image the superficial vasculature; the deep
    # capillary-only plexuses are essentially invisible to them (OCTA sees
    # them instead), so the HRF comparison rasterizes layer 0 only
    fundus = edges[edges.layer_id == 0]
    raster = metrics.rasterize_network(fundus, geometry.bounds, radii=fundus.radius.values)
    image = metrics.image_metrics(raster)
    stats = image_stats(image, references)
    # Bifurcation geometry is judged on the superficial tree, like the raster.
    # Angles are measured in 3D on the tree; the plexus is nearly planar, so
    # they agree with the fundus (x-y) projection the literature reports
    angles = metrics.bifurcation_angles(pop[pop.layer_id == 0])
    superficial = pop[pop.layer_id == 0]

    def perfused(vessels, vessel_type=None) -> float:
        return metrics.perfused_fraction(
            vessels, geometry.semi_axes, *geometry.perfusion, vessel_type=vessel_type
        )

    stats.update(
        {
            "bifurcation_angle_median": (
                float(np.median(angles)) if len(angles) else float("nan")
            ),
            "bifurcation_obtuse_share": (
                float((angles > 100).mean()) if len(angles) else float("nan")
            ),
            "artery_vein_caliber_ratio": metrics.arcade_caliber_ratio(
                pop, geometry.disc_center
            ),
            # Scored perfusion needs both supply and drainage; the any-vessel
            # and per-tree fractions are the growth diagnostics behind it
            "perfused_fraction": metrics.paired_perfused_fraction(
                pop, geometry.semi_axes, *geometry.perfusion
            ),
            "colonized_fraction": perfused(pop),
            "superficial_colonized_fraction": perfused(superficial),
            "arterial_supply_fraction": perfused(pop, metrics.VESSEL_TYPE_ARTERY),
            "venous_drainage_fraction": perfused(pop, metrics.VESSEL_TYPE_VEIN),
            "fov_fraction": image["fov_fraction"],
            "n_frozen": int(pop.frozen.sum()),
            "graph_cycles": metrics.graph_cycles(pop),
        }
    )
    return stats


def apply_overrides(spec: dict, overrides: dict) -> dict:
    candidate = copy.deepcopy(spec)
    for (section, key), value in overrides.items():
        candidate["configuration"][section][key] = value
    return candidate


def evaluate_spec(spec: dict, steps: int, references: dict, workdir: Path, tag: str) -> dict:
    """Run one candidate spec to completion and score it."""
    sim = simulation.build_from_spec(spec, workdir / f"candidate_{tag}.yaml")
    geometry = simulation.get_geometry(sim)
    simulation.run_steps(sim, steps)
    pop = simulation.get_network(sim)
    edges = simulation.tree_edges(pop)
    stats = scoring_stats(pop, edges, geometry, references)
    return {"stats": stats, "scores": calibration_score(stats)}


def coordinate_descent(evaluate, space: dict, base_overrides: dict, budget: int) -> dict:
    """Greedy one-knob-at-a-time search; ``evaluate(overrides) -> total score``.

    Returns {"best": overrides, "best_score": float, "evaluations": int}.
    The evaluate callable is expected to cache/log as it sees fit.
    """
    current = dict(base_overrides)
    best_score = evaluate(current)
    evaluations = 1
    improved = True
    while improved and evaluations < budget:
        improved = False
        for knob, candidates in space.items():
            for value in candidates:
                if evaluations >= budget:
                    break
                if current.get(knob) == value:
                    continue
                trial = {**current, knob: value}
                score = evaluate(trial)
                evaluations += 1
                if score < best_score:
                    best_score = score
                    current = trial
                    improved = True
    return {"best": current, "best_score": best_score, "evaluations": evaluations}


def combine_seed_scores(per_seed: dict) -> dict:
    """Mean of each score component across seeds — the multi-seed objective.

    Averaging the totals (rather than taking the best or median seed) makes
    a config that collapses on any seed lose to one that is merely mediocre
    everywhere, which is the robustness the single-seed fits kept missing.
    """
    combined: dict[str, float] = {}
    for scores in per_seed.values():
        for name, value in scores.items():
            combined[name] = combined.get(name, 0.0) + value / len(per_seed)
    return combined


@click.command()
@click.argument("model_spec", type=click.Path(exists=True))
@click.option("--budget", default=24, show_default=True, help="Max simulation runs.")
@click.option("--steps", default=800, show_default=True)
@click.option("--seed", default=None, type=int, help="Override the spec's random seed.")
@click.option(
    "--seeds",
    default=None,
    help="Comma-separated random seeds; the objective becomes the MEAN score "
    "across them (one simulation per seed per evaluation), so the fit "
    "cannot win by overfitting a lucky growth trajectory.",
)
@click.option(
    "--log",
    "log_path",
    default="calibration_log.json",
    show_default=True,
    type=click.Path(),
)
@click.option(
    "--workdir",
    default=".calibration_work",
    show_default=True,
    type=click.Path(),
    help="Scratch directory for candidate specs (not for committing).",
)
@click.option(
    "--derive-targets",
    is_flag=True,
    help="Print the HRF-derived TARGETS (mean/sd across masks) under the "
    "current measurement conventions and exit, without simulating.",
)
def main(
    model_spec: str,
    budget: int,
    steps: int,
    seed,
    seeds,
    log_path: str,
    workdir: str,
    derive_targets: bool,
):
    """Fit MODEL_SPEC's knobs against the HRF-derived calibration targets."""
    click.echo("Computing HRF reference statistics (cached download)...")
    references = hrf_references()
    if derive_targets:
        for name, spec in derive_hrf_targets(references).items():
            current = TARGETS[name]
            click.echo(
                f"  {name:28s} target {spec['target']:.4f}  scale {spec['scale']:.4f}"
                f"   (TARGETS: {current['target']:.4f} / {current['scale']:.4f})"
            )
        return

    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    with open(model_spec) as f:
        base_spec = yaml.safe_load(f)
    if seed is not None:
        base_spec = simulation.with_seed(base_spec, seed)
    seed_list = [int(s) for s in seeds.split(",")] if seeds else None

    log: list[dict] = []
    cache: dict[tuple, float] = {}
    counter = itertools.count()

    def evaluate(overrides: dict) -> float:
        key = tuple(sorted(overrides.items()))
        if key in cache:
            return cache[key]
        tag = f"{next(counter):03d}"
        start = time.time()
        candidate = apply_overrides(base_spec, overrides)
        if seed_list is None:
            result = evaluate_spec(candidate, steps, references, work, tag)
            entry = dict(result)
            total = result["scores"]["total"]
        else:
            per_seed = {}
            for s in seed_list:
                per_seed[s] = evaluate_spec(
                    simulation.with_seed(candidate, s),
                    steps,
                    references,
                    work,
                    f"{tag}_seed{s}",
                )
            combined = combine_seed_scores({s: r["scores"] for s, r in per_seed.items()})
            entry = {
                "scores": combined,
                "per_seed": {str(s): r for s, r in per_seed.items()},
            }
            total = combined["total"]
        cache[key] = total
        log.append(
            {
                "tag": tag,
                "overrides": {f"{s}.{k}": v for (s, k), v in overrides.items()},
                "minutes": round((time.time() - start) / 60, 1),
                **entry,
            }
        )
        with open(log_path, "w") as f:
            json.dump(log, f, indent=2)
        click.echo(f"[{tag}] total={total:.2f} {log[-1]['overrides']}")
        return total

    # Seed the search at the spec's own settings so the incumbent value of
    # each knob is recognized and never re-evaluated as a "new" candidate
    base_overrides = {
        (section, key): base_spec["configuration"][section][key]
        for section, key in SEARCH_SPACE
        if key in base_spec["configuration"].get(section, {})
    }
    outcome = coordinate_descent(evaluate, SEARCH_SPACE, base_overrides, budget)
    click.echo(f"Best score {outcome['best_score']:.2f} after {outcome['evaluations']} runs:")
    for (section, key), value in outcome["best"].items():
        click.echo(f"  {section}.{key} = {value}")
    click.echo(f"Full log in {log_path}")


if __name__ == "__main__":
    main()

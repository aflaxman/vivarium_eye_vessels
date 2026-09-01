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
# HRF-derived targets use the across-mask mean and sd; clinical targets use
# literature values with a judgment scale; "one_sided" components only score
# deviations in the bad direction.
TARGETS = {
    "skeleton_density": {"target": 0.0321, "scale": 0.0029},
    "area_density": {"target": 0.1193, "scale": 0.0103},
    "fractal_dimension": {"target": 1.3536, "scale": 0.0242},
    "branch_tortuosity_median": {"target": 1.000, "scale": 0.02},
    "ks_log_length": {"target": 0.0, "scale": 0.05, "one_sided": "above"},
    "capillary_share": {"target": 0.0647, "scale": 0.05, "one_sided": "above"},
    # HRF branches are 34.5% wide (>4 px); shares sum to 1, so scoring the
    # capillary and wide shares also pins the mid (2-4 px) share
    "wide_share": {"target": 0.345, "scale": 0.08},
    # Real arcades run straight: the 90th-percentile tortuosity of wide
    # (>4 px) branches is 1.11 across the HRF masks (sd 0.017). One-sided
    # with a judgment scale: only meandering wide vessels are penalized
    "wide_tortuosity_q90": {"target": 1.11, "scale": 0.05, "one_sided": "above"},
    # Comb-like side branching: real arcades carry a branch point every
    # ~23 px of wide (>4 px) skeleton (HRF across-mask mean 22.71, sd 1.77,
    # counting connected junction clusters once)
    "wide_junction_spacing_px": {"target": 22.71, "scale": 1.77},
    # Fundus-visible (superficial-tree) bifurcation geometry: healthy
    # arteriolar branch angles are unimodal around ~75-84 degrees; obtuse
    # (>100 degree) junctions are rare. Judgment scales from the literature
    "bifurcation_angle_median": {"target": 77.0, "scale": 5.0},
    "bifurcation_obtuse_share": {"target": 0.05, "scale": 0.05, "one_sided": "above"},
    "artery_vein_caliber_ratio": {"target": 0.67, "scale": 0.05},
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
    ("frozen_repulsion", "interaction_radius"): [0.08, 0.1, 0.15],
    ("path_freezer", "radius_taper"): [0.994, 0.996, 0.998],
    ("flow_remodeler", "shear_threshold_fraction"): [0.35, 0.5, 0.65],
    ("flow_remodeler", "adaptation_rate"): [0.05, 0.10, 0.15],
    ("flow_remodeler", "adaptation_deadband"): [1.0, 2.0, 4.0],
    ("flow_remodeler", "max_radius"): [0.012, 0.016, 0.02],
    ("flow_remodeler", "max_adapted_radius"): [0.005, 0.006, 0.008],
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


def pooled_hrf_lengths() -> np.ndarray:
    lengths: list[float] = []
    for path in reference_data.fetch_hrf_masks():
        binary = metrics.binarize_mask(reference_data.load_mask(path))
        lengths.extend(metrics.image_metrics(binary)["branch_length_px"])
    return np.asarray(lengths, dtype=float)


def scoring_stats(pop, edges, bounds, semi_axes, real_lengths: np.ndarray) -> dict:
    """The scored summary statistics for one finished simulation."""
    from vivarium_eye_vessels.vnv.compare import stratify_by_diameter

    # Fundus photographs image the superficial vasculature; the deep
    # capillary-only plexuses are essentially invisible to them (OCTA sees
    # them instead), so the HRF comparison rasterizes layer 0 only
    fundus = edges[edges.layer_id == 0]
    raster = metrics.rasterize_network(fundus, bounds, radii=fundus.radius.values)
    image = metrics.image_metrics(raster)
    lengths = np.asarray(image["branch_length_px"], dtype=float)
    strata = stratify_by_diameter(lengths, image["branch_diameter_px"])
    tortuosity = np.asarray(image["branch_tortuosity"], dtype=float)
    diameter = np.asarray(image["branch_diameter_px"], dtype=float)
    wide_tortuosity = tortuosity[diameter > 4.0]
    # Bifurcation geometry is judged on the superficial tree, like the raster
    angles = metrics.bifurcation_angles(pop[pop.layer_id == 0])
    arteries = pop[(pop.vessel_type == 1) & (pop.depth == 0) & (pop.radius > 0)]
    veins = pop[(pop.vessel_type == 2) & (pop.depth == 0) & (pop.radius > 0)]
    return {
        "skeleton_density": image["skeleton_density"],
        "area_density": image["area_density"],
        "fractal_dimension": image["fractal_dimension"],
        "branch_tortuosity_median": (
            float(np.median(image["branch_tortuosity"]))
            if image["branch_tortuosity"]
            else float("nan")
        ),
        "ks_log_length": (
            float(ks_2samp(np.log10(lengths), np.log10(real_lengths)).statistic)
            if len(lengths)
            else float("nan")
        ),
        "capillary_share": len(strata["diameter_le_2px"]) / max(len(lengths), 1),
        "wide_share": len(strata["diameter_gt_4px"]) / max(len(lengths), 1),
        "wide_tortuosity_q90": (
            float(np.quantile(wide_tortuosity, 0.9)) if len(wide_tortuosity) else float("nan")
        ),
        "wide_junction_spacing_px": image["wide_junction_spacing_px"],
        "bifurcation_angle_median": (
            float(np.median(angles)) if len(angles) else float("nan")
        ),
        "bifurcation_obtuse_share": (
            float((angles > 100).mean()) if len(angles) else float("nan")
        ),
        "artery_vein_caliber_ratio": (
            float(arteries.radius.mean() / veins.radius.mean())
            if len(arteries) and len(veins)
            else float("nan")
        ),
        "perfused_fraction": metrics.perfused_fraction(pop, semi_axes, 0.1, 0.15),
        "n_frozen": int(pop.frozen.sum()),
        "graph_cycles": metrics.graph_cycles(pop),
    }


def apply_overrides(spec: dict, overrides: dict) -> dict:
    candidate = copy.deepcopy(spec)
    for (section, key), value in overrides.items():
        candidate["configuration"][section][key] = value
    return candidate


def evaluate_spec(
    spec: dict, steps: int, real_lengths: np.ndarray, workdir: Path, tag: str
) -> dict:
    """Run one candidate spec to completion and score it."""
    spec_path = workdir / f"candidate_{tag}.yaml"
    with open(spec_path, "w") as f:
        yaml.safe_dump(spec, f)
    sim = simulation.build_headless_simulation(spec_path)
    bounds = simulation.get_ellipsoid_bounds(sim)
    semi_axes = simulation.get_ellipsoid_semi_axes(sim)
    simulation.run_steps(sim, steps)
    pop = simulation.get_network(sim)
    edges = simulation.tree_edges(pop)
    stats = scoring_stats(pop, edges, bounds, semi_axes, real_lengths)
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
def main(model_spec: str, budget: int, steps: int, seed, seeds, log_path: str, workdir: str):
    """Fit MODEL_SPEC's knobs against the HRF-derived calibration targets."""
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    with open(model_spec) as f:
        base_spec = yaml.safe_load(f)
    if seed is not None:
        base_spec["configuration"]["randomness"]["random_seed"] = seed
    seed_list = [int(s) for s in seeds.split(",")] if seeds else None

    click.echo("Computing HRF reference statistics (cached download)...")
    real_lengths = pooled_hrf_lengths()

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
            result = evaluate_spec(candidate, steps, real_lengths, work, tag)
            entry = dict(result)
            total = result["scores"]["total"]
        else:
            per_seed = {}
            for s in seed_list:
                candidate["configuration"]["randomness"]["random_seed"] = s
                per_seed[s] = evaluate_spec(
                    candidate, steps, real_lengths, work, f"{tag}_seed{s}"
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

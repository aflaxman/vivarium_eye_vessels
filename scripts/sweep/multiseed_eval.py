"""Evaluate one set of spec overrides on one or more seeds; write scores and stats.

Usage: python scripts/sweep/multiseed_eval.py OUTDIR NAME '{"section.key": value}' SEED[,SEED...]

Runs the model spec with the dotted-key overrides applied, 800 steps per
seed, and scores each run with :mod:`vivarium_eye_vessels.vnv.calibrate`.
Prints a one-line summary per seed, the combined mean, and finally the
``MULTI-DONE`` marker that sweep_jobs.sh treats as completion; writes
``OUTDIR/NAME_multiseed.json`` with every seed's scores and statistics.
"""

import json
import sys
import time
from pathlib import Path

import yaml

from vivarium_eye_vessels.vnv import calibrate

SPEC = (
    Path(__file__).resolve().parents[2]
    / "src/vivarium_eye_vessels/model_specifications/model_spec.yaml"
)


def summary(name: str, seed: int, result: dict, minutes: float) -> str:
    scores, stats = result["scores"], result["stats"]
    return (
        f"{name} seed{seed}: total={scores['total']:.1f} n={stats['n_frozen']}"
        f" perf={stats['perfused_fraction']:.2f} art={stats['arterial_supply_fraction']:.2f}"
        f" ven={stats['venous_drainage_fraction']:.2f}"
        f" skel={stats['skeleton_density'] * 100:.2f}% area={stats['area_density'] * 100:.2f}%"
        f" FD={stats['fractal_dimension']:.3f} ksCal={stats['ks_caliber_profile']:.3f}"
        f" spacing={stats['wide_junction_spacing_px']:.1f}"
        f" align={stats.get('arcade_radial_alignment', float('nan')):.2f}"
        f" reach={stats.get('arcade_reach_px', float('nan')):.0f}"
        f" thick={stats.get('thick_share', float('nan')):.3f}"
        f" avr={stats['artery_vein_caliber_ratio']:.2f} ({minutes:.1f}m)"
    )


def main() -> None:
    out = Path(sys.argv[1])
    name = sys.argv[2]
    overrides = json.loads(sys.argv[3])
    seeds = [int(seed) for seed in sys.argv[4].split(",")]
    out.mkdir(parents=True, exist_ok=True)

    references = calibrate.hrf_references()
    with open(SPEC) as f:
        base_spec = yaml.safe_load(f)

    per_seed: dict[int, dict] = {}
    for seed in seeds:
        start = time.time()
        spec = calibrate.apply_dotted_overrides(base_spec, overrides)
        spec["configuration"]["randomness"]["random_seed"] = seed
        result = calibrate.evaluate_spec(spec, 800, references, out, f"{name}_s{seed}")
        per_seed[seed] = result
        print(summary(name, seed, result, (time.time() - start) / 60), flush=True)
        with open(out / f"{name}_multiseed.json", "w") as f:
            json.dump(per_seed, f, indent=2, default=str)

    combined = calibrate.combine_seed_scores(
        {seed: result["scores"] for seed, result in per_seed.items()}
    )
    print(f"{name} MEAN total={combined['total']:.1f}", flush=True)
    print("MULTI-DONE", flush=True)


if __name__ == "__main__":
    main()

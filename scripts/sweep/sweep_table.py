"""Per-config mean score decomposition and mean statistics across seeds.

Usage: python scripts/sweep/sweep_table.py [NAME ...]
Reads $SWEEP_DIR/NAME_sSEED_multiseed.json (default sweep_out/); with no
names, every config found there.
"""

import glob
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

SWEEP_DIR = Path(os.environ.get("SWEEP_DIR", "sweep_out"))
STATS = [
    "perfused_fraction",
    "arterial_supply_fraction",
    "venous_drainage_fraction",
    "colonized_fraction",
    "skeleton_density",
    "area_density",
    "fractal_dimension",
    "capillary_share",
    "wide_share",
    "wide_junction_spacing_px",
    "ks_caliber_profile",
    "artery_vein_caliber_ratio",
    "bifurcation_obtuse_share",
    "branch_tortuosity_mean",
    "arcade_radial_alignment",
    "arcade_reach_px",
    "thick_share",
    "n_frozen",
]


def main() -> None:
    names = sys.argv[1:] or sorted(
        {
            Path(path).name.rsplit("_s", 1)[0]
            for path in glob.glob(str(SWEEP_DIR / "*_multiseed.json"))
        }
    )
    per = defaultdict(dict)
    for name in names:
        for path in glob.glob(str(SWEEP_DIR / f"{name}_s[0-9]*_multiseed.json")):
            for seed, result in json.load(open(path)).items():
                per[name][int(seed)] = result
    names = [name for name in names if per[name]]
    if not names:
        print("no results found in", SWEEP_DIR)
        return
    keys = [k for k in next(iter(per[names[0]].values()))["scores"] if k != "total"]
    print(f"{'mean score':26s}" + "".join(f"{name[:13]:>14}" for name in names))
    for key in keys:
        print(
            f"{key:26s}"
            + "".join(
                f"{np.mean([r['scores'][key] for r in per[name].values()]):14.1f}"
                for name in names
            )
        )
    print(
        f"{'TOTAL mean':26s}"
        + "".join(
            f"{np.mean([r['scores']['total'] for r in per[name].values()]):14.1f}"
            for name in names
        )
    )

    def per_seed(name: str) -> str:
        return "/".join(f"{per[name][s]['scores']['total']:.0f}" for s in sorted(per[name]))

    print(f"{'per seed':26s}" + "".join(f"{per_seed(name):>14}" for name in names))
    print()
    for key in STATS:
        print(
            f"{key:26s}"
            + "".join(
                f"{np.mean([r['stats'].get(key, np.nan) for r in per[name].values()]):14.4f}"
                for name in names
            )
        )


if __name__ == "__main__":
    main()

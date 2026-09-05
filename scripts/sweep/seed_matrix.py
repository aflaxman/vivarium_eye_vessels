"""Per-seed matrix of total score and perfusion vitals for named sweep configs.

Usage: python scripts/sweep/seed_matrix.py NAME [NAME ...]
Reads $SWEEP_DIR/NAME_sSEED_multiseed.json (default sweep_out/).
"""

import glob
import json
import os
import sys
from pathlib import Path

SWEEP_DIR = Path(os.environ.get("SWEEP_DIR", "sweep_out"))


def load(names: list[str]) -> dict:
    data = {}
    for name in names:
        for path in glob.glob(str(SWEEP_DIR / f"{name}_s[0-9]*_multiseed.json")):
            for seed, result in json.load(open(path)).items():
                data[(name, int(seed))] = result
    return data


def main() -> None:
    names = sys.argv[1:]
    data = load(names)
    seeds = sorted({seed for _, seed in data})
    print(f"{'config':16s}" + "".join(f"{seed:>22}" for seed in seeds) + f"{'mean':>10}")
    for name in names:
        cells, totals = [], []
        for seed in seeds:
            result = data.get((name, seed))
            if result is None:
                cells.append(f"{'-':>22}")
                continue
            stats, scores = result["stats"], result["scores"]
            totals.append(scores["total"])
            cells.append(
                f"{scores['total']:6.0f} p{stats['perfused_fraction']:.2f}"
                f" a{stats['arterial_supply_fraction']:.2f} v{stats['venous_drainage_fraction']:.2f}"
            )
        mean = sum(totals) / len(totals) if totals else float("nan")
        print(f"{name:16s}" + "".join(cells) + f"{mean:10.1f} (n={len(totals)})")


if __name__ == "__main__":
    main()

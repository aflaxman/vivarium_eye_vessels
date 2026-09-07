"""Score run_sim_dump.py pickles against the current calibrate.TARGETS.

Usage: python scripts/sweep/score_dump.py DUMP.pkl [DUMP.pkl ...]

Re-scores finished runs after a change to the targets or the estimators
without re-simulating; prints the terms worth half a point or more.
"""

import pickle
import sys

from vivarium_eye_vessels.vnv import calibrate


def main() -> None:
    references = calibrate.hrf_references()
    for path in sys.argv[1:]:
        with open(path, "rb") as f:
            dump = pickle.load(f)
        stats = calibrate.scoring_stats(
            dump["pop"], dump["edges"], dump["geometry"], references
        )
        scores = calibrate.calibration_score(stats)
        print(path, f"total={scores['total']:.1f}")
        for key, value in sorted(scores.items(), key=lambda item: -item[1]):
            if key != "total" and value >= 0.5:
                print(f"   {key:28s} score={value:6.1f}  value={stats.get(key):.4g}")


if __name__ == "__main__":
    main()

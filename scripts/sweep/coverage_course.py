"""Time course of per-tree coverage, paired perfusion, live tips and tip deaths.

Usage: python scripts/sweep/coverage_course.py ['{"section.key": value}'] [SEED]

Prints one row every 50 steps: wave radius, arterial / venous / paired
coverage of the demand lattice (all sites and the far periphery), live
tips per tree with their median distance from the disc and how many are
at the branching depth limit, and the new dead ends per tree since the
previous row. This is the diagnostic that located the paired-perfusion
deficit (roadmap eighteenth pass).
"""

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.spatial import cKDTree

from vivarium_eye_vessels.components.boundaries import generate_demand_sites
from vivarium_eye_vessels.components.particles import (
    VESSEL_TYPE_ARTERY,
    VESSEL_TYPE_VEIN,
)
from vivarium_eye_vessels.vnv import calibrate, simulation

SPEC = (
    Path(__file__).resolve().parents[2]
    / "src/vivarium_eye_vessels/model_specifications/model_spec.yaml"
)
SWEEP_DIR = Path(os.environ.get("SWEEP_DIR", "sweep_out"))
COLUMNS = [
    "x",
    "y",
    "z",
    "frozen",
    "vessel_type",
    "path_id",
    "depth",
    "anastomosis_id",
    "parent_id",
]


def main() -> None:
    overrides = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    seed = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != "-" else None
    with open(SPEC) as f:
        spec = calibrate.apply_dotted_overrides(yaml.safe_load(f), overrides)
    if seed is not None:
        spec = simulation.with_seed(spec, seed)
    SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    sim = simulation.build_from_spec(spec, SWEEP_DIR / "coverage_course.yaml")
    geometry = simulation.get_geometry(sim)
    spacing, perfusion_radius = geometry.perfusion
    disc = np.asarray(geometry.disc_center)
    sites = generate_demand_sites(np.asarray(geometry.semi_axes, float), spacing)
    far = np.hypot(sites[:, 0] - disc[0], sites[:, 1] - disc[1]) > 2.0
    max_depth = int(spec["configuration"]["path_splitter"]["max_depth"])
    wave = sim.get_component("developmental_wave")
    rows: list[dict] = []
    previous_dead = {VESSEL_TYPE_ARTERY: set(), VESSEL_TYPE_VEIN: set()}

    def reach(frozen: pd.DataFrame, vessel_type: int) -> np.ndarray:
        sub = frozen[frozen.vessel_type == vessel_type][["x", "y", "z"]].to_numpy()
        if len(sub) == 0:
            return np.zeros(len(sites), bool)
        return cKDTree(sub).query(sites, k=1)[0] <= perfusion_radius

    def snapshot(step: int) -> None:
        pop = sim.get_population(COLUMNS)
        frozen = pop[pop.frozen & (pop.vessel_type > 0)]
        arterial, venous = reach(frozen, VESSEL_TYPE_ARTERY), reach(frozen, VESSEL_TYPE_VEIN)
        children = set(pop.parent_id[pop.parent_id >= 0])
        row = {
            "step": step,
            "radius": round(wave.radius, 2) if wave is not None else np.nan,
            "art": arterial.mean(),
            "ven": venous.mean(),
            "paired": (arterial & venous).mean(),
            "art_far": arterial[far].mean(),
            "ven_far": venous[far].mean(),
        }
        for vessel_type, label in ((VESSEL_TYPE_ARTERY, "A"), (VESSEL_TYPE_VEIN, "V")):
            sub = pop[pop.vessel_type == vessel_type]
            tips = sub[~sub.frozen & (sub.path_id >= 0)]
            row[f"tips_{label}"] = len(tips)
            row[f"tips_{label}_maxdepth"] = int((tips.depth >= max_depth).sum())
            row[f"tipdist_{label}"] = (
                float(np.median(np.hypot(tips.x - disc[0], tips.y - disc[1])))
                if len(tips)
                else np.nan
            )
            dead = sub[sub.frozen & ~sub.index.isin(children) & (sub.anastomosis_id < 0)]
            new_dead = set(dead.index) - previous_dead[vessel_type]
            previous_dead[vessel_type] = set(dead.index)
            fresh = dead.loc[list(new_dead)]
            row[f"newdead_{label}"] = len(fresh)
            row[f"newdead_{label}_far"] = (
                int((np.hypot(fresh.x - disc[0], fresh.y - disc[1]) > 2.0).sum())
                if len(fresh)
                else 0
            )
        rows.append(row)

    snapshot(0)
    simulation.run_steps(sim, 800, callback=snapshot, every=50)
    pd.set_option("display.width", 300)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()

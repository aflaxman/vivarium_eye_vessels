"""Run the model spec with overrides and pickle the finished network for diagnostics.

Usage: python scripts/sweep/run_sim_dump.py OUT.pkl ['{"section.key": value}'] [SEED] [STEPS]

The pickle carries the particle table, the tree edges, the superficial
raster, the containment bounds and the V&V geometry, so score_dump.py and
ad-hoc diagnostics can work on a finished run without re-simulating.
"""

import json
import pickle
import sys
import time
from pathlib import Path

import yaml

from vivarium_eye_vessels.vnv import calibrate, metrics, simulation

SPEC = (
    Path(__file__).resolve().parents[2]
    / "src/vivarium_eye_vessels/model_specifications/model_spec.yaml"
)


def main() -> None:
    out = Path(sys.argv[1])
    overrides = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    seed = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != "-" else None
    steps = int(sys.argv[4]) if len(sys.argv) > 4 else 800

    with open(SPEC) as f:
        spec = calibrate.apply_dotted_overrides(yaml.safe_load(f), overrides)
    if seed is not None:
        spec = simulation.with_seed(spec, seed)

    start = time.time()
    sim = simulation.build_from_spec(spec, out.with_suffix(".yaml"))
    geometry = simulation.get_geometry(sim)
    simulation.run_steps(sim, steps)
    pop = sim.get_population(
        [
            "x",
            "y",
            "z",
            "vx",
            "vy",
            "vz",
            "frozen",
            "freeze_time",
            "unfreeze_time",
            "depth",
            "parent_id",
            "path_id",
            "radius",
            "vessel_type",
            "anastomosis_id",
            "layer_id",
        ]
    )
    edges = simulation.tree_edges(pop)
    fundus = edges[edges.layer_id == 0]
    raster = metrics.rasterize_network(fundus, geometry.bounds, radii=fundus.radius.values)
    remodeler = sim.get_component("flow_remodeler")
    with open(out, "wb") as f:
        pickle.dump(
            {
                "pop": pop,
                "edges": edges,
                "raster": raster,
                "geometry": geometry,
                "bounds": geometry.bounds,
                "semi_axes": geometry.semi_axes,
                "total_pruned": (
                    int(remodeler.total_pruned) if remodeler is not None else None
                ),
                "overrides": overrides,
                "seed": seed,
                "seconds": time.time() - start,
            },
            f,
        )
    print(
        f"dumped {out} in {time.time() - start:.0f} s; n_frozen {int(pop.frozen.sum())}",
        f"pruned {int(remodeler.total_pruned) if remodeler is not None else None}",
    )


if __name__ == "__main__":
    main()

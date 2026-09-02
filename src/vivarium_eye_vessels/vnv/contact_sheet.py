"""Seed-reliability contact sheet: fresh-seed simulations vs HRF masks.

One healthy retina looks much like another, so a realistic model must
produce a usable network on *every* seed, not just the calibration seeds.
This sheet runs the model on held-out seeds, renders each superficial
network beside expert-labeled HRF masks in the same style, and stamps each
simulation panel with its reliability vitals (perfused fraction, skeleton
density, frozen count). The JSON alongside records the same numbers so
seed reliability is tracked across model versions like every other metric.

Usage::

    vnv_contact_sheet src/vivarium_eye_vessels/model_specifications/model_spec.yaml
"""

import json
import time
from pathlib import Path

import click
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from vivarium_eye_vessels.vnv import metrics, reference_data, simulation

DEFAULT_SEEDS = "11,202,909,4242"
DEFAULT_MASKS = "03_h.tif,06_h.tif,09_h.tif,14_h.tif"
PERFUSION_TARGET = 0.95  # a seed counts as reliable at or above this


def run_seed(spec: dict, seed: int, steps: int, workdir: Path) -> dict:
    """Run one seed and return its superficial raster and reliability vitals."""
    candidate = yaml.safe_load(yaml.safe_dump(spec))
    candidate["configuration"]["randomness"]["random_seed"] = seed
    spec_path = workdir / f"contact_seed{seed}.yaml"
    with open(spec_path, "w") as f:
        yaml.safe_dump(candidate, f)
    sim = simulation.build_headless_simulation(spec_path)
    bounds = simulation.get_ellipsoid_bounds(sim)
    semi_axes = simulation.get_ellipsoid_semi_axes(sim)
    simulation.run_steps(sim, steps)
    pop = simulation.get_network(sim)
    edges = simulation.tree_edges(pop)
    fundus = edges[edges.layer_id == 0]
    raster = metrics.rasterize_network(fundus, bounds, radii=fundus.radius.values)
    return {
        "seed": seed,
        "raster": raster,
        "perfused_fraction": metrics.perfused_fraction(pop, semi_axes, 0.1, 0.15),
        "skeleton_density": metrics.image_metrics(raster)["skeleton_density"],
        "n_frozen": int(pop.frozen.sum()),
    }


@click.command()
@click.argument("model_spec", type=click.Path(exists=True))
@click.option("--seeds", default=DEFAULT_SEEDS, show_default=True)
@click.option("--masks", default=DEFAULT_MASKS, show_default=True)
@click.option("--steps", default=800, show_default=True)
@click.option(
    "--output-dir",
    default="docs/vnv",
    show_default=True,
    type=click.Path(),
)
def main(model_spec: str, seeds: str, masks: str, steps: int, output_dir: str):
    """Render held-out-seed simulations beside HRF masks for MODEL_SPEC."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with open(model_spec) as f:
        spec = yaml.safe_load(f)
    seed_list = [int(s) for s in seeds.split(",")]
    mask_names = [m.strip() for m in masks.split(",")]

    runs = []
    for seed in seed_list:
        start = time.time()
        runs.append(run_seed(spec, seed, steps, output))
        run = runs[-1]
        click.echo(
            f"seed {seed}: perfused {run['perfused_fraction']:.2f}, "
            f"skeleton {run['skeleton_density']*100:.2f}%, "
            f"n_frozen {run['n_frozen']} "
            f"({(time.time() - start) / 60:.1f} min)"
        )
        (output / f"contact_seed{seed}.yaml").unlink(missing_ok=True)

    columns = max(len(runs), len(mask_names))
    fig, axes = plt.subplots(2, columns, figsize=(5.5 * columns, 10))
    axes = np.atleast_2d(axes)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, run in zip(axes[0], runs):
        ax.imshow(~run["raster"], cmap="gray", interpolation="nearest")
        reliable = run["perfused_fraction"] >= PERFUSION_TARGET
        ax.set_title(
            f"sim seed {run['seed']} — perfused {run['perfused_fraction']:.0%}, "
            f"skeleton {run['skeleton_density']*100:.1f}%"
            + ("" if reliable else "  [STALLED]"),
            fontsize=11,
            color="black" if reliable else "firebrick",
        )
    available = {p.name: p for p in reference_data.fetch_hrf_masks()}
    for ax, name in zip(axes[1], mask_names):
        binary = metrics.binarize_mask(reference_data.load_mask(available[name]))
        ax.imshow(~binary, cmap="gray", interpolation="nearest")
        ax.set_title(f"HRF {name}", fontsize=11)

    reliability = float(
        np.mean([run["perfused_fraction"] >= PERFUSION_TARGET for run in runs])
    )
    fig.suptitle(
        f"Held-out seeds (top) vs HRF expert masks (bottom) — "
        f"{reliability:.0%} of seeds reach {PERFUSION_TARGET:.0%} perfusion",
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    sheet_path = output / "contact_sheet.png"
    fig.savefig(sheet_path, dpi=100)

    record = {
        "seeds": [
            {key: value for key, value in run.items() if key != "raster"} for run in runs
        ],
        "perfusion_target": PERFUSION_TARGET,
        "seed_reliability": reliability,
        "steps": steps,
        "masks": mask_names,
    }
    with open(output / "contact_sheet.json", "w") as f:
        json.dump(record, f, indent=2)
    click.echo(f"Seed reliability: {reliability:.0%}. Wrote {sheet_path}")


if __name__ == "__main__":
    main()

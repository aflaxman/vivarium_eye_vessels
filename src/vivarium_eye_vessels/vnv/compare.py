"""Diagnostic comparison of simulated networks against real vessel masks.

Usage::

    vnv_compare src/vivarium_eye_vessels/model_specifications/model_spec.yaml

Runs the model specification headless, computes network metrics on the result,
computes the same image-based metrics on expert-labeled vessel masks from the
public HRF dataset (downloaded to a local cache on first use), and writes:

- ``comparison.png``: side-by-side visual and distributional diagnostics
- ``metrics.json``: every metric plus run parameters, for quantitative
  tracking across model versions

By default the outputs overwrite ``docs/vnv/`` in place, so committing them
makes the before/after of a model change visible as image and JSON diffs in
the pull request.
"""

import datetime
import json
from pathlib import Path

import click
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from skimage.morphology import dilation, disk, skeletonize

from vivarium_eye_vessels.vnv import metrics, reference_data, simulation

SIM_COLOR = "#2a78d6"  # categorical slot 1
REAL_COLOR = "#eb6834"  # categorical slot 2
INK = "#333333"
MUTED = "#767676"
GRID = "#e6e6e6"
LITERATURE_ANGLE_RANGE = (60, 90)  # healthy retinal bifurcations (Masters 2004)


def style_axis(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def overlay_histogram(ax, sim_values, real_values, bins, xlabel: str) -> None:
    """Overlaid density histograms: simulation vs. real, step outlines."""
    common = dict(bins=bins, density=True, histtype="step", linewidth=1.8)
    if len(real_values):
        ax.hist(real_values, color=REAL_COLOR, label="HRF (real)", **common)
    if len(sim_values):
        ax.hist(sim_values, color=SIM_COLOR, label="Simulation", **common)
    ax.set_xlabel(xlabel, color=INK, fontsize=9)
    ax.set_ylabel("density", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    style_axis(ax)


def run_comparison(model_spec: str, output_dir: Path, steps: int) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Simulate ---
    sim = simulation.build_headless_simulation(model_spec)
    bounds = simulation.get_ellipsoid_bounds(sim)
    semi_axes = simulation.get_ellipsoid_semi_axes(sim)
    simulation.run_steps(sim, steps)
    pop = simulation.get_network(sim)
    edges = simulation.tree_edges(pop)

    has_calibers = "radius" in edges.columns and bool((edges.radius > 0).any())
    sim_raster = metrics.rasterize_network(
        edges, bounds, radii=edges.radius.values if has_calibers else None
    )
    sim_image_metrics = metrics.image_metrics(sim_raster)
    sim_angles = metrics.bifurcation_angles(pop)
    sim_tortuosity_paths = metrics.path_tortuosity(pop)
    sim_junction_exponents = metrics.junction_exponents(pop)
    sim_perfused_fraction = metrics.perfused_fraction(
        pop, semi_axes, site_spacing=0.1, perfusion_radius=0.15
    )
    sim_arterial_supply = metrics.perfused_fraction(
        pop, semi_axes, site_spacing=0.1, perfusion_radius=0.15, vessel_type=1
    )
    sim_venous_drainage = metrics.perfused_fraction(
        pop, semi_axes, site_spacing=0.1, perfusion_radius=0.15, vessel_type=2
    )
    arteries = pop[(pop.vessel_type == 1) & (pop.path_id >= 0) & (pop.radius > 0)]
    veins = pop[(pop.vessel_type == 2) & (pop.path_id >= 0) & (pop.radius > 0)]
    # Clinical AVR is measured on the major arcades near the disc, so compare
    # the depth-0 trunks rather than averaging over capillary-scale tails
    arcade_arteries = arteries[arteries.depth == 0]
    arcade_veins = veins[veins.depth == 0]
    if len(arcade_arteries) and len(arcade_veins):
        sim_avr = float(arcade_arteries.radius.mean() / arcade_veins.radius.mean())
    else:
        sim_avr = float("nan")

    # --- Real reference data ---
    mask_paths = reference_data.fetch_hrf_masks()
    real_per_mask = []
    real_lengths: list[float] = []
    real_tortuosity: list[float] = []
    example_binary = None
    for path in mask_paths:
        binary = metrics.binarize_mask(reference_data.load_mask(path))
        if example_binary is None:
            example_binary = binary
        m = metrics.image_metrics(binary)
        m["file"] = path.name
        real_lengths.extend(m.pop("branch_length_px"))
        real_tortuosity.extend(m.pop("branch_tortuosity"))
        real_per_mask.append(m)

    real_fd = np.array([m["fractal_dimension"] for m in real_per_mask])
    real_density = np.array([m["skeleton_density"] for m in real_per_mask])
    real_area_density = np.array([m["area_density"] for m in real_per_mask])

    # --- Figure ---
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.patch.set_facecolor("white")

    def displayable(thin_image: np.ndarray) -> np.ndarray:
        # 1px skeletons disappear when matplotlib downsamples; thicken for display only
        return dilation(thin_image, disk(2))

    ax = axes[0, 0]
    sim_display = sim_raster if has_calibers else displayable(sim_raster)
    sim_panel_title = (
        "Simulation: network with calibers (x–y)"
        if has_calibers
        else "Simulation: rasterized network (x–y)"
    )
    ax.imshow(~sim_display, cmap="gray", interpolation="nearest")
    ax.set_title(sim_panel_title, color=INK, fontsize=10)
    for spine in ax.spines.values():
        spine.set_color(SIM_COLOR)
        spine.set_linewidth(2)
    ax.set_xticks([])
    ax.set_yticks([])

    ax = axes[0, 1]
    ax.imshow(~example_binary, cmap="gray", interpolation="nearest")
    ax.set_title(
        f"HRF healthy eye ({mask_paths[0].name}): expert-labeled vessels",
        color=INK,
        fontsize=10,
    )
    for spine in ax.spines.values():
        spine.set_color(REAL_COLOR)
        spine.set_linewidth(2)
    ax.set_xticks([])
    ax.set_yticks([])

    ax = axes[0, 2]
    ax.imshow(~displayable(skeletonize(example_binary)), cmap="gray", interpolation="nearest")
    ax.set_title("HRF skeleton (as measured)", color=INK, fontsize=10)
    for spine in ax.spines.values():
        spine.set_color(REAL_COLOR)
        spine.set_linewidth(2)
    ax.set_xticks([])
    ax.set_yticks([])

    length_bins = np.geomspace(
        metrics.MIN_BRANCH_PIXELS,
        max(
            max(real_lengths, default=100),
            max(sim_image_metrics["branch_length_px"], default=100),
        ),
        30,
    )
    ax = axes[1, 0]
    overlay_histogram(
        ax,
        sim_image_metrics["branch_length_px"],
        real_lengths,
        length_bins,
        "skeleton branch length (px)",
    )
    ax.set_xscale("log")
    ax.set_title("Segment length distribution", color=INK, fontsize=10)

    ax = axes[1, 1]
    tortuosity_bins = np.linspace(1.0, 2.0, 30)
    overlay_histogram(
        ax,
        np.clip(sim_image_metrics["branch_tortuosity"], 1, 2),
        np.clip(real_tortuosity, 1, 2),
        tortuosity_bins,
        "branch tortuosity (arc / chord)",
    )
    ax.set_title("Tortuosity distribution", color=INK, fontsize=10)

    ax = axes[1, 2]
    if len(sim_angles):
        ax.hist(
            sim_angles,
            bins=np.linspace(0, 180, 25),
            density=True,
            histtype="step",
            linewidth=1.8,
            color=SIM_COLOR,
            label="Simulation",
        )
    ax.axvspan(
        *LITERATURE_ANGLE_RANGE,
        color=REAL_COLOR,
        alpha=0.12,
        label="literature range (healthy retina)",
    )
    ax.set_xlabel("bifurcation angle (degrees)", color=INK, fontsize=9)
    ax.set_ylabel("density", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    ax.set_title("Bifurcation angles (tree-based)", color=INK, fontsize=10)
    style_axis(ax)
    if len(sim_junction_exponents):
        ax.text(
            0.97,
            0.55,
            f"junction exponent:\nmedian k = {np.median(sim_junction_exponents):.2f}"
            f" (n={len(sim_junction_exponents)})\nMurray target k = 3",
            transform=ax.transAxes,
            ha="right",
            va="top",
            color=INK,
            fontsize=8,
        )

    headline = (
        f"Fractal dimension (skeleton): sim {sim_image_metrics['fractal_dimension']:.2f}  vs  "
        f"HRF {real_fd.mean():.2f} ± {real_fd.std():.2f}      "
        f"Skeleton density: sim {sim_image_metrics['skeleton_density']*100:.2f}%  vs  "
        f"HRF {real_density.mean()*100:.2f}% ± {real_density.std()*100:.2f}%"
    )
    area_line = (
        f"Vessel area density: sim {sim_image_metrics['area_density']*100:.2f}%  vs  "
        f"HRF {real_area_density.mean()*100:.2f}% ± {real_area_density.std()*100:.2f}%      "
        f"Perfused tissue: sim {sim_perfused_fraction*100:.1f}%      "
        f"A:V caliber ratio: sim {sim_avr:.2f} (clinical ~0.67)"
    )
    fig.suptitle(
        "Vessel network diagnostics: simulation vs. HRF public dataset",
        color=INK,
        fontsize=13,
        y=0.995,
    )
    fig.text(0.5, 0.955, headline, ha="center", color=INK, fontsize=10)
    fig.text(0.5, 0.93, area_line, ha="center", color=INK, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.915))
    figure_path = output_dir / "comparison.png"
    fig.savefig(figure_path, dpi=110, facecolor="white")
    plt.close(fig)

    # --- Metrics JSON ---
    results = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "model_spec": str(model_spec),
        "steps": steps,
        "raster_size": metrics.RASTER_SIZE,
        "has_calibers": has_calibers,
        "simulation": {
            "n_particles": int(len(pop)),
            "n_frozen": int(pop.frozen.sum()),
            "n_segments": int(len(edges)),
            "fractal_dimension": sim_image_metrics["fractal_dimension"],
            "skeleton_density": sim_image_metrics["skeleton_density"],
            "area_density": sim_image_metrics["area_density"],
            "branch_length_px": metrics.summarize(sim_image_metrics["branch_length_px"]),
            "branch_tortuosity": metrics.summarize(sim_image_metrics["branch_tortuosity"]),
            "path_tortuosity": metrics.summarize(sim_tortuosity_paths),
            "bifurcation_angle_deg": metrics.summarize(sim_angles),
            "junction_exponent": metrics.summarize(sim_junction_exponents),
            "perfused_fraction": sim_perfused_fraction,
            "arterial_supply_fraction": sim_arterial_supply,
            "venous_drainage_fraction": sim_venous_drainage,
            "n_artery_segments": int(len(arteries)),
            "n_vein_segments": int(len(veins)),
            "artery_vein_caliber_ratio": sim_avr,
        },
        "real_hrf": {
            "n_masks": len(real_per_mask),
            "fractal_dimension": {
                "mean": float(real_fd.mean()),
                "std": float(real_fd.std()),
            },
            "skeleton_density": {
                "mean": float(real_density.mean()),
                "std": float(real_density.std()),
            },
            "area_density": {
                "mean": float(real_area_density.mean()),
                "std": float(real_area_density.std()),
            },
            "branch_length_px": metrics.summarize(real_lengths),
            "branch_tortuosity": metrics.summarize(real_tortuosity),
            "per_mask": real_per_mask,
        },
    }
    metrics_path = output_dir / "metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(results, f, indent=2)

    click.echo(f"Wrote {figure_path} and {metrics_path}")
    return results


@click.command()
@click.argument("model_spec", type=click.Path(exists=True))
@click.option("-o", "--output-dir", default="docs/vnv", show_default=True, type=click.Path())
@click.option("--steps", default=800, show_default=True, help="Total simulation steps.")
def main(model_spec: str, output_dir: str, steps: int) -> None:
    """Compare MODEL_SPEC's vessel network to real HRF vessel masks."""
    run_comparison(model_spec, Path(output_dir), steps)


if __name__ == "__main__":
    main()

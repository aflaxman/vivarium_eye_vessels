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
import time
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
DIAMETER_BIN_EDGES = (2.0, 4.0)  # px: capillary (<=2), mid (2-4], wide (>4)
DIAMETER_BIN_TITLES = {
    "diameter_le_2px": "diameter ≤ 2 px (capillary)",
    "diameter_2_4px": "diameter 2–4 px",
    "diameter_gt_4px": "diameter > 4 px (arcades)",
}


def stratify_by_diameter(lengths, diameters) -> dict[str, np.ndarray]:
    """Branch lengths split into capillary / mid / wide diameter strata.

    Diameters come from the distance transform along the skeleton, so the
    same stratification applies to the sim raster and to real masks that
    carry no explicit calibers. The capillary bin is closed at 2 px because
    that is the transform's floor: a single-pixel line measures exactly 2.
    """
    lengths = np.asarray(lengths, dtype=float)
    diameters = np.asarray(diameters, dtype=float)
    lo, hi = DIAMETER_BIN_EDGES
    return {
        "diameter_le_2px": lengths[diameters <= lo],
        "diameter_2_4px": lengths[(diameters > lo) & (diameters <= hi)],
        "diameter_gt_4px": lengths[diameters > hi],
    }


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


LAYER_NAMES = ["superficial", "intermediate", "deep"]
LAYER_COLORS = ["#2a78d6", "#eb6834", "#8a63c9"]


def layer_name(layer: int) -> str:
    return LAYER_NAMES[layer] if layer < len(LAYER_NAMES) else f"layer {layer}"


def render_plexus_figure(pop, edges, bounds, layer_z, output_path: Path) -> None:
    """En-face slab per plexus plus an x-z cross-section — an OCTA-style view."""
    from matplotlib.collections import LineCollection

    layers = sorted(int(layer) for layer in edges.layer_id.unique() if layer >= 0)
    n_panels = len(layers) + 1
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5.4))
    fig.patch.set_facecolor("white")

    for column, layer in enumerate(layers):
        ax = axes[column]
        in_layer = edges[edges.layer_id == layer]
        raster = metrics.rasterize_network(in_layer, bounds, radii=in_layer.radius.values)
        ax.imshow(~raster, cmap="gray", interpolation="nearest")
        ax.set_title(
            f"{layer_name(layer)} plexus (z = {layer_z[layer]:+.2f}) — "
            f"{len(in_layer)} segments",
            color=INK,
            fontsize=10,
        )
        for spine in ax.spines.values():
            spine.set_color(LAYER_COLORS[layer % len(LAYER_COLORS)])
            spine.set_linewidth(2)
        ax.set_xticks([])
        ax.set_yticks([])

    # Cross-section: the stratification is invisible in the fundus view
    ax = axes[-1]
    segments = np.stack(
        [edges[["x0", "z0"]].to_numpy(float), edges[["x1", "z1"]].to_numpy(float)],
        axis=1,
    )
    colors = [
        LAYER_COLORS[int(layer) % len(LAYER_COLORS)] if layer >= 0 else MUTED
        for layer in edges.layer_id
    ]
    ax.add_collection(LineCollection(segments, colors=colors, linewidths=0.7, alpha=0.6))
    for layer in layers:
        ax.axhline(
            layer_z[layer],
            color=LAYER_COLORS[layer % len(LAYER_COLORS)],
            linewidth=0.8,
            linestyle="--",
            alpha=0.7,
        )
    a, _ = bounds
    ax.set_xlim(-a * 1.05, a * 1.05)
    z_extent = max(abs(min(layer_z)), abs(max(layer_z))) * 2.5
    ax.set_ylim(-z_extent, z_extent)
    ax.set_xlabel("x", color=INK, fontsize=9)
    ax.set_ylabel("z", color=INK, fontsize=9)
    ax.set_title("Cross-section (x–z): the layers themselves", color=INK, fontsize=10)
    style_axis(ax)

    fig.suptitle(
        "Stratified plexuses: en-face slabs and cross-section",
        color=INK,
        fontsize=13,
        y=0.99,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output_path, dpi=110, facecolor="white")
    plt.close(fig)


def plexus_metrics(pop, edges, layer_z) -> dict:
    """Per-plexus composition and stratification quality."""
    vessels = pop[(pop.layer_id >= 0) & (pop.radius > 0)]
    per_layer = []
    for layer in sorted(int(layer) for layer in vessels.layer_id.unique()):
        in_layer = vessels[vessels.layer_id == layer]
        per_layer.append(
            {
                "layer": layer,
                "name": layer_name(layer),
                "z_plane": float(layer_z[layer]),
                "n_segments": int(len(in_layer)),
                "max_radius": float(in_layer.radius.max()),
                "median_abs_z_error": float((in_layer.z - layer_z[layer]).abs().median()),
            }
        )
    on_path = pop[pop.layer_id >= 0]
    children = on_path[on_path.parent_id.isin(on_path.index)]
    n_diving = int(
        (
            children.layer_id.to_numpy()
            != on_path.layer_id.loc[children.parent_id].to_numpy()
        ).sum()
    )
    return {
        "layer_z": [float(z) for z in layer_z],
        "per_layer": per_layer,
        "n_diving_vessels": n_diving,
    }


def run_comparison(model_spec: str, output_dir: Path, steps: int) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Simulate ---
    setup_start = time.perf_counter()
    sim = simulation.build_headless_simulation(model_spec)
    setup_seconds = time.perf_counter() - setup_start
    bounds = simulation.get_ellipsoid_bounds(sim)
    semi_axes = simulation.get_ellipsoid_semi_axes(sim)
    run_start = time.perf_counter()
    simulation.run_steps(sim, steps)
    simulation_seconds = time.perf_counter() - run_start
    pop = simulation.get_network(sim)
    edges = simulation.tree_edges(pop)

    has_calibers = "radius" in edges.columns and bool((edges.radius > 0).any())
    # Fundus photographs image the superficial vasculature; the deep
    # capillary-only plexuses are essentially invisible to them (OCTA sees
    # them instead — docs/vnv/plexus.png), so everything compared against
    # HRF uses the superficial (layer 0) projection only
    fundus_edges = edges[edges.layer_id == 0] if "layer_id" in edges.columns else edges
    sim_raster = metrics.rasterize_network(
        fundus_edges, bounds, radii=fundus_edges.radius.values if has_calibers else None
    )
    sim_image_metrics = metrics.image_metrics(sim_raster)
    # Angles and Murray exponents are compared against literature values for
    # fundus-visible (arteriolar) junctions, so measure them on the
    # superficial tree — the deep capillary plexuses form polygonal meshes
    # whose T-shaped junctions would drown the arteriolar geometry
    superficial_pop = pop[pop.layer_id == 0]
    sim_angles = metrics.bifurcation_angles(superficial_pop)
    sim_angles_all_layers = metrics.bifurcation_angles(pop)
    sim_tortuosity_paths = metrics.path_tortuosity(pop)
    sim_tree_segment_lengths = metrics.tree_segment_lengths(pop)
    sim_n_anastomoses = int((pop.anastomosis_id >= 0).sum())
    sim_graph_cycles = metrics.graph_cycles(pop)
    remodeler = sim.get_component("flow_remodeler")
    sim_n_pruned = int(remodeler.total_pruned) if remodeler is not None else 0
    plexus = sim.get_component("plexus_layers")
    sim_plexus = None
    if plexus is not None and "layer_id" in pop.columns:
        sim_plexus = plexus_metrics(pop, edges, plexus.layer_z)
        render_plexus_figure(pop, edges, bounds, plexus.layer_z, output_dir / "plexus.png")
    sim_shear = (
        remodeler.solve_network(sim.get_population(remodeler.required_attributes))
        if remodeler is not None
        else None
    )
    sim_junction_exponents = metrics.junction_exponents(superficial_pop)
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
    real_diameters: list[float] = []
    real_pixel_diameters: list[np.ndarray] = []
    example_binary = None
    for path in mask_paths:
        binary = metrics.binarize_mask(reference_data.load_mask(path))
        real_pixel_diameters.append(metrics.skeleton_pixel_diameters(binary))
        if example_binary is None:
            example_binary = binary
        m = metrics.image_metrics(binary)
        m["file"] = path.name
        real_lengths.extend(m.pop("branch_length_px"))
        real_tortuosity.extend(m.pop("branch_tortuosity"))
        real_diameters.extend(m.pop("branch_diameter_px"))
        real_per_mask.append(m)

    sim_strata = stratify_by_diameter(
        sim_image_metrics["branch_length_px"], sim_image_metrics["branch_diameter_px"]
    )
    real_strata = stratify_by_diameter(real_lengths, real_diameters)
    real_pixel_diameters = np.concatenate(real_pixel_diameters)
    sim_pixel_diameters = metrics.skeleton_pixel_diameters(sim_raster)

    # Calibration score: squared z-like deviation from each validation target
    from scipy.stats import ks_2samp

    from vivarium_eye_vessels.vnv import calibrate

    sim_lengths = np.asarray(sim_image_metrics["branch_length_px"], dtype=float)
    sim_branch_tortuosity = np.asarray(sim_image_metrics["branch_tortuosity"], dtype=float)
    sim_branch_diameter = np.asarray(sim_image_metrics["branch_diameter_px"], dtype=float)
    sim_wide_tortuosity = sim_branch_tortuosity[sim_branch_diameter > 4.0]
    calibration_stats = {
        "skeleton_density": sim_image_metrics["skeleton_density"],
        "area_density": sim_image_metrics["area_density"],
        "fractal_dimension": sim_image_metrics["fractal_dimension"],
        "branch_tortuosity_median": (
            float(np.median(sim_image_metrics["branch_tortuosity"]))
            if sim_image_metrics["branch_tortuosity"]
            else float("nan")
        ),
        "ks_log_length": (
            float(
                ks_2samp(np.log10(sim_lengths), np.log10(np.asarray(real_lengths))).statistic
            )
            if len(sim_lengths)
            else float("nan")
        ),
        "capillary_share": len(sim_strata["diameter_le_2px"]) / max(len(sim_lengths), 1),
        "wide_share": len(sim_strata["diameter_gt_4px"]) / max(len(sim_lengths), 1),
        "wide_tortuosity_q90": (
            float(np.quantile(sim_wide_tortuosity, 0.9))
            if len(sim_wide_tortuosity)
            else float("nan")
        ),
        "wide_junction_spacing_px": sim_image_metrics["wide_junction_spacing_px"],
        "ks_caliber_profile": (
            float(ks_2samp(sim_pixel_diameters, real_pixel_diameters).statistic)
            if len(sim_pixel_diameters)
            else float("nan")
        ),
        "bifurcation_angle_median": (
            float(np.median(sim_angles)) if len(sim_angles) else float("nan")
        ),
        "bifurcation_obtuse_share": (
            float((sim_angles > 100).mean()) if len(sim_angles) else float("nan")
        ),
        "artery_vein_caliber_ratio": sim_avr,
        "perfused_fraction": sim_perfused_fraction,
    }
    calibration_scores = calibrate.calibration_score(calibration_stats)

    real_fd = np.array([m["fractal_dimension"] for m in real_per_mask])
    real_density = np.array([m["skeleton_density"] for m in real_per_mask])
    real_area_density = np.array([m["area_density"] for m in real_per_mask])

    # --- Figure ---
    fig, axes = plt.subplots(4, 3, figsize=(15, 18))
    fig.patch.set_facecolor("white")

    def displayable(thin_image: np.ndarray) -> np.ndarray:
        # 1px skeletons disappear when matplotlib downsamples; thicken for display only
        return dilation(thin_image, disk(2))

    ax = axes[0, 0]
    sim_display = sim_raster if has_calibers else displayable(sim_raster)
    sim_panel_title = (
        "Simulation: superficial network with calibers (x–y)"
        if has_calibers
        else "Simulation: rasterized superficial network (x–y)"
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
    ax.set_title("Bifurcation angles (superficial, tree-based)", color=INK, fontsize=10)
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

    # Row 3: segment lengths stratified by local vessel diameter, measured the
    # same way on both images (distance transform along the skeleton)
    sim_n_branches = max(len(sim_lengths), 1)
    real_n_branches = max(len(real_lengths), 1)
    for column, key in enumerate(DIAMETER_BIN_TITLES):
        ax = axes[2, column]
        overlay_histogram(
            ax,
            sim_strata[key],
            real_strata[key],
            length_bins,
            "skeleton branch length (px)",
        )
        ax.set_xscale("log")
        ax.set_title(
            f"{DIAMETER_BIN_TITLES[key]} — sim {len(sim_strata[key])/sim_n_branches:.0%},"
            f" HRF {len(real_strata[key])/real_n_branches:.0%} of branches",
            color=INK,
            fontsize=10,
        )

    # Row 4: the diameter composition itself — the distribution and the
    # share of branches (and of total skeleton length) per diameter stratum
    ax = axes[3, 0]
    # Length-weighted caliber profile: one sample per skeleton pixel at its
    # local (EDT) diameter, so the histogram is skeleton length by width.
    # EDT diameters are quantized (2, 2.83, 4, 4.47, ...), so 1 px bins
    overlay_histogram(
        ax,
        sim_pixel_diameters,
        real_pixel_diameters,
        np.arange(1.5, 13.6, 1.0),
        "local vessel diameter at skeleton px",
    )
    ax.set_title("Caliber profile (skeleton length by width)", color=INK, fontsize=10)

    def composition_bars(ax, sim_shares, real_shares, ylabel: str, title: str) -> None:
        labels = ["≤ 2 px", "2–4 px", "> 4 px"]
        positions = np.arange(len(labels))
        width = 0.36
        ax.bar(positions - width / 2, sim_shares, width, color=SIM_COLOR, label="Simulation")
        ax.bar(
            positions + width / 2, real_shares, width, color=REAL_COLOR, label="HRF (real)"
        )
        for xs, shares in (
            (positions - width / 2, sim_shares),
            (positions + width / 2, real_shares),
        ):
            for x, share in zip(xs, shares):
                ax.text(x, share + 0.01, f"{share:.0%}", ha="center", color=INK, fontsize=8)
        ax.set_xticks(positions, labels)
        ax.set_ylim(0, 1.12)
        ax.set_ylabel(ylabel, color=INK, fontsize=9)
        ax.legend(frameon=False, fontsize=8, labelcolor=INK)
        ax.set_title(title, color=INK, fontsize=10)
        style_axis(ax)

    stratum_keys = list(DIAMETER_BIN_TITLES)
    composition_bars(
        axes[3, 1],
        [len(sim_strata[key]) / sim_n_branches for key in stratum_keys],
        [len(real_strata[key]) / real_n_branches for key in stratum_keys],
        "share of branches",
        "Diameter composition (by branch count)",
    )
    sim_total_length = max(float(np.sum(sim_lengths)), 1.0)
    real_total_length = max(float(np.sum(real_lengths)), 1.0)
    composition_bars(
        axes[3, 2],
        [float(np.sum(sim_strata[key])) / sim_total_length for key in stratum_keys],
        [float(np.sum(real_strata[key])) / real_total_length for key in stratum_keys],
        "share of skeleton length",
        "Diameter composition (by skeleton length)",
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
        f"A:V caliber ratio: sim {sim_avr:.2f} (clinical ~0.67)      "
        f"Loops: {sim_graph_cycles}      Pruned: {sim_n_pruned}      "
        f"Calibration score: {calibration_scores['total']:.1f}"
    )
    fig.suptitle(
        "Vessel network diagnostics: simulation vs. HRF public dataset",
        color=INK,
        fontsize=13,
        y=0.995,
    )
    fig.text(0.5, 0.976, headline, ha="center", color=INK, fontsize=10)
    fig.text(0.5, 0.963, area_line, ha="center", color=INK, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
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
        # Wall-clock runtime on the machine that generated this file; useful
        # for tracking the trend across model versions, not as an absolute
        "runtime": {
            "setup_seconds": round(setup_seconds, 2),
            "simulation_seconds": round(simulation_seconds, 2),
            "steps_per_second": round(steps / simulation_seconds, 2),
        },
        "simulation": {
            "n_particles": int(len(pop)),
            "n_frozen": int(pop.frozen.sum()),
            "n_segments": int(len(edges)),
            "fractal_dimension": sim_image_metrics["fractal_dimension"],
            "skeleton_density": sim_image_metrics["skeleton_density"],
            "area_density": sim_image_metrics["area_density"],
            "branch_length_px": metrics.summarize(sim_image_metrics["branch_length_px"]),
            "branch_diameter_px": metrics.summarize(sim_image_metrics["branch_diameter_px"]),
            "branch_length_by_diameter": {
                key: {
                    **metrics.summarize(values),
                    "share": len(values) / max(len(sim_image_metrics["branch_length_px"]), 1),
                    "length_share": float(np.sum(values))
                    / max(float(np.sum(sim_lengths)), 1.0),
                }
                for key, values in sim_strata.items()
            },
            "tree_segment_length": metrics.summarize(sim_tree_segment_lengths),
            "branch_tortuosity": metrics.summarize(sim_image_metrics["branch_tortuosity"]),
            "path_tortuosity": metrics.summarize(sim_tortuosity_paths),
            "pixel_diameter_px": metrics.summarize(sim_pixel_diameters),
            "bifurcation_angle_deg": metrics.summarize(sim_angles),
            "bifurcation_angle_deg_all_layers": metrics.summarize(sim_angles_all_layers),
            "junction_exponent": metrics.summarize(sim_junction_exponents),
            "perfused_fraction": sim_perfused_fraction,
            "arterial_supply_fraction": sim_arterial_supply,
            "venous_drainage_fraction": sim_venous_drainage,
            "n_artery_segments": int(len(arteries)),
            "n_vein_segments": int(len(veins)),
            "artery_vein_caliber_ratio": sim_avr,
            "n_anastomoses": sim_n_anastomoses,
            "graph_cycles": sim_graph_cycles,
            "n_pruned": sim_n_pruned,
            "plexus_layers": sim_plexus,
            "calibration": {
                "stats": calibration_stats,
                "scores": calibration_scores,
            },
            "wall_shear": (
                metrics.summarize(sim_shear.shear[sim_shear.shear > 0])
                if sim_shear is not None
                else metrics.summarize([])
            ),
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
            "branch_diameter_px": metrics.summarize(real_diameters),
            "branch_length_by_diameter": {
                key: {
                    **metrics.summarize(values),
                    "share": len(values) / max(len(real_lengths), 1),
                    "length_share": float(np.sum(values))
                    / max(float(np.sum(np.asarray(real_lengths, dtype=float))), 1.0),
                }
                for key, values in real_strata.items()
            },
            "branch_tortuosity": metrics.summarize(real_tortuosity),
            "pixel_diameter_px": metrics.summarize(real_pixel_diameters),
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

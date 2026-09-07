"""Quantitative vessel network metrics for validation against real data.

Two families of metrics:

- **Tree-based** (simulation only): computed directly from the particle
  table's parent/child structure — bifurcation angles and per-path tortuosity.
- **Image-based** (simulation and real masks): computed on a binary skeleton
  image so that the caliber-less simulated centerlines and the expert-labeled
  vessel masks are measured the same way — box-counting fractal dimension,
  skeleton density, per-branch segment lengths and tortuosity.

Spatial scale in the simulation is arbitrary, so image-based metrics are
computed at a common raster resolution and reported in pixels or as fractions
of the imaged region. They support before/after comparison across model
versions, not absolute physiological claims.

Measurement conventions, applied identically to both sources so that the
comparison is apples to apples:

- **Binarization is a majority vote.** A pixel is vessel when more than half
  of it is covered. HRF masks are block-averaged down to the common
  resolution and thresholded at 0.5; simulated segments are drawn by the
  exact pixel-center rule that is the limit of drawing at infinite
  resolution and majority-downsampling. Vessels narrower than half a pixel
  vanish from both, as they do from a fundus photograph.
- **Densities are per imaged pixel**, the convex hull of the vessel pixels,
  not per frame pixel: a fundus camera's circular field and the model's
  elliptical territory both leave empty corners.
- **The skeleton is pruned of spurs** shorter than ``MIN_BRANCH_PIXELS``
  before anything is counted, so ragged edges do not manufacture branch
  points, and branch length is the arc length of the pixel chain (diagonal
  steps count sqrt 2), so a straight 45-degree line has tortuosity 1.
- **Box counting uses a fixed range of box sizes** (2..128 px) rather than
  one derived from the frame, so padding a skeleton into a larger canvas
  does not change its fractal dimension.
- **Arcade geometry is read relative to an estimated disc**, for both
  sources: the point the wide vessels' tangent lines converge on
  (:func:`disc_center`). A simulated raster is windowed to the reference
  image's extent first (:func:`fundus_window`), since how far the arcades
  reach depends on how far the field extends from the disc.
- **The macula is the largest vessel-free disk near the fovea**
  (:func:`clear_radius_px`), at fundus scale on the skeleton and at OCTA
  scale on the thresholded en-face signal (:func:`faz_metrics`), because an
  inscribed disk cannot leak through a gap the way a region can. One
  simulation unit is ``MM_PER_UNIT`` (4.5 mm), the pin that lets OCTA scans
  and the model be drawn at the same micrometers per pixel.
"""

from typing import Any

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.optimize import brentq
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import cKDTree
from skimage.draw import line as draw_line
from skimage.measure import perimeter
from skimage.morphology import convex_hull_image, skeletonize

from vivarium_eye_vessels.components.particles import (
    VESSEL_TYPE_ARTERY,
    VESSEL_TYPE_VEIN,
)

RASTER_SIZE = 1024
MIN_BRANCH_PIXELS = 5
# Distance band from the disc center (simulation units) in which the
# artery:vein caliber ratio is read off the arcades -- the analog of the
# clinical measurement zone around the disc margin
AVR_ZONE = (0.1, 0.5)
BOX_SIZES = 2 ** np.arange(1, 8)  # 2..128 px; both sources span at least 512 px
NEIGHBOR_KERNEL = np.ones((3, 3), dtype=int)
NEIGHBOR_KERNEL[1, 1] = 0


##################
# Raster helpers #
##################


def rasterize_network(
    edges: pd.DataFrame,
    bounds: tuple[float, float],
    size: int = RASTER_SIZE,
    radii: np.ndarray | None = None,
) -> np.ndarray:
    """Rasterize vessel segments into a binary image (x-y projection).

    Parameters
    ----------
    edges
        Segment endpoints as produced by
        :func:`~vivarium_eye_vessels.vnv.simulation.tree_edges`.
    bounds
        The (a, b) semi-axes of the containment region; the image spans
        [-a, a] x [-b, b].
    size
        Output image size in pixels (longest dimension).
    radii
        Optional per-segment vessel radii (in simulation units, aligned with
        ``edges``). When given, a pixel is vessel when its center lies within
        a segment's radius — the majority-vote rule the HRF masks are
        binarized by (see the module docstring). Without radii, segments are
        drawn one pixel wide.
    """
    a, b = bounds
    aspect = b / a
    width = size if a >= b else max(int(size * a / b), 1)
    height = max(int(width * aspect), 1)

    image = np.zeros((height, width), dtype=bool)
    if edges.empty:
        return image

    # Endpoints in continuous pixel coordinates (pixel centers at integers)
    scale_x, scale_y = (width - 1) / (2 * a), (height - 1) / (2 * b)
    c0, c1 = (edges.x0.values + a) * scale_x, (edges.x1.values + a) * scale_x
    r0, r1 = (edges.y0.values + b) * scale_y, (edges.y1.values + b) * scale_y

    if radii is None:
        rows = np.clip(np.round([r0, r1]).astype(int), 0, height - 1)
        cols = np.clip(np.round([c0, c1]).astype(int), 0, width - 1)
        for i in range(len(edges)):
            rr, cc = draw_line(rows[0, i], cols[0, i], rows[1, i], cols[1, i])
            image[rr, cc] = True
        return image

    # A tube of radius r covers more than half of a pixel centered at
    # distance d when d < r (r >= 0.5 px); a thinner tube only where
    # d <= 0.5 - r, and one narrower than half a pixel nowhere
    radii_px = np.asarray(radii, dtype=float) * scale_x
    reach = np.where(radii_px >= 0.5, radii_px, np.where(radii_px > 0.25, 0.5 - radii_px, 0))
    for i in np.nonzero(reach > 0)[0]:
        rmin, rmax = max(int(np.floor(min(r0[i], r1[i]) - reach[i])), 0), min(
            int(np.ceil(max(r0[i], r1[i]) + reach[i])), height - 1
        )
        cmin, cmax = max(int(np.floor(min(c0[i], c1[i]) - reach[i])), 0), min(
            int(np.ceil(max(c0[i], c1[i]) + reach[i])), width - 1
        )
        if rmax < rmin or cmax < cmin:
            continue
        rows = np.arange(rmin, rmax + 1)[:, None]
        cols = np.arange(cmin, cmax + 1)[None, :]
        dr, dc = r1[i] - r0[i], c1[i] - c0[i]
        length2 = dr * dr + dc * dc
        along = (
            np.clip(((rows - r0[i]) * dr + (cols - c0[i]) * dc) / length2, 0.0, 1.0)
            if length2 > 0
            else 0.0
        )
        distance2 = (rows - r0[i] - along * dr) ** 2 + (cols - c0[i] - along * dc) ** 2
        image[rmin : rmax + 1, cmin : cmax + 1] |= distance2 < reach[i] ** 2
    return image


def fundus_window(
    raster: np.ndarray, reference_shape: tuple[int, int]
) -> tuple[np.ndarray, tuple[int, int, int, int]]:
    """Central crop of ``raster`` with the pixel extent of a reference image.

    Returns the crop and its (row0, col0, rows, cols) placement. Because the
    simulation is rasterized at the same pixels-per-caliber scale as the HRF
    masks, a window the size of the HRF working image shows the two at the
    same magnification and aspect, and statistics that depend on how far the
    field extends from the disc (arcade reach) are read over the same extent.
    """
    rows = min(reference_shape[0], raster.shape[0])
    cols = min(reference_shape[1], raster.shape[1])
    row0 = (raster.shape[0] - rows) // 2
    col0 = (raster.shape[1] - cols) // 2
    return raster[row0 : row0 + rows, col0 : col0 + cols], (row0, col0, rows, cols)


def rasterize_window(
    edges: pd.DataFrame,
    center: tuple[float, float],
    width_units: float,
    size_px: int,
    radii: np.ndarray | None = None,
) -> np.ndarray:
    """Rasterize a square window of the network, ``width_units`` wide around ``center``.

    Same drawing rule as :func:`rasterize_network`; segments that cannot
    touch the window are dropped first so the border is not drawn on.
    """
    half = width_units / 2
    shifted = edges.assign(
        x0=edges.x0 - center[0],
        x1=edges.x1 - center[0],
        y0=edges.y0 - center[1],
        y1=edges.y1 - center[1],
    )
    margin = 0.0 if radii is None else float(np.max(radii, initial=0.0))
    inside = (
        (np.minimum(shifted.x0, shifted.x1) <= half + margin)
        & (np.maximum(shifted.x0, shifted.x1) >= -half - margin)
        & (np.minimum(shifted.y0, shifted.y1) <= half + margin)
        & (np.maximum(shifted.y0, shifted.y1) >= -half - margin)
    ).to_numpy()
    kept_radii = None if radii is None else np.asarray(radii)[inside]
    return rasterize_network(shifted[inside], (half, half), size=size_px, radii=kept_radii)


####################
# Physical scale   #
####################

# The simulation's length unit is pinned to the retina through the fundus
# raster: HRF calibers were matched in pixels at RASTER_SIZE over the 4-unit
# field, and the model's disc-to-fovea distance (1.05 units) is the anatomical
# 4.76 mm, so one unit is 4.5 mm. The root vein caliber (0.034 units) then
# reads 153 um, a central retinal vein equivalent, and the HRF working image
# (876 px) spans 15 mm, a 45-degree fundus: the pin is consistent to ~5%
MM_PER_UNIT = 4.5
FUNDUS_PX_PER_UNIT = RASTER_SIZE / 4.0  # shared by the HRF working raster and the sim raster
FUNDUS_MM_PER_PX = MM_PER_UNIT / FUNDUS_PX_PER_UNIT
FUNDUS_WINDOW_SHAPE = (584, 876)  # the HRF working image: 2336 x 3504 block-averaged by 4
# OCTA: the ROSE-1 en-face scans are 3 x 3 mm at 304 px, centered on the fovea
OCTA_SCAN_MM = 3.0
OCTA_SIZE_PX = 304
OCTA_MM_PER_PX = OCTA_SCAN_MM / OCTA_SIZE_PX
MACULA_SEARCH_MM = 1.5  # a macula-centered image has its fovea within this of the center
FAZ_SMOOTH_MM = 0.04  # en-face signal is smoothed at about one capillary spacing
FAZ_THRESHOLD_FRACTION = 0.5  # vascular at or above this fraction of the scan's mean signal
FAZ_SEARCH_MM = 0.3  # the FAZ's clear disk is centered within this of the image center
FAZ_SEED_MM = 0.1  # the zone must overlap a disk this size at the image center
# The capillary class: CapillaryBed sprouts are 0.0009 units (8 um) wide, strictly
# below the 0.001 adaptation floor of the arteriole tree. Junction and path
# statistics compared against fundus references skip them, as a fundus does
CAPILLARY_RADIUS_UNITS = 0.001


def vascular_signal(signal: np.ndarray, mm_per_px: float) -> np.ndarray:
    """Where an en-face signal carries vessels: at or above half its smoothed mean.

    Reads a real OCTA angiogram (flow signal, capillaries bright) and a
    simulated network drawn as 0/1 pixels the same way: the signal is
    smoothed at ``FAZ_SMOOTH_MM`` (about one capillary spacing) and pixels
    at or above ``FAZ_THRESHOLD_FRACTION`` of the smoothed signal's mean
    over the image count as vascular. The mean, not the median, so a sparse
    network drawn on an empty background (most of the model's macular
    window is background at capillary resolution) still has a threshold
    above zero; the whole scan, not its center, so a zone that swallows the
    center still measures.
    """
    smooth = ndimage.gaussian_filter(
        np.asarray(signal, dtype=float), FAZ_SMOOTH_MM / mm_per_px
    )
    return smooth >= FAZ_THRESHOLD_FRACTION * smooth.mean()


def clear_radius_px(vascular: np.ndarray, search_radius_px: float) -> float:
    """Radius (px) of the largest vessel-free disk centered near the image center.

    The distance from the nearest vascular pixel, maximized over centers
    within ``search_radius_px`` of the image center, which is where a
    macula-centered image and the simulation's windows both put the fovea.
    An inscribed disk cannot leak through a gap the way a region can, so
    the same statistic serves the fundus (the vessel skeleton as the
    vascular mask: the macula a photograph shows clear) and OCTA (the
    thresholded flow signal: the FAZ). NaN when nothing is vascular.
    """
    vascular = np.asarray(vascular, dtype=bool)
    if not vascular.any():
        return float("nan")
    clearance = ndimage.distance_transform_edt(~vascular)
    rows, cols = np.indices(vascular.shape)
    near = (
        np.hypot(rows - vascular.shape[0] / 2, cols - vascular.shape[1] / 2)
        <= search_radius_px
    )
    return float(clearance[near].max())


def avascular_zone(signal: np.ndarray, mm_per_px: float) -> np.ndarray:
    """The foveal avascular zone of an en-face signal as a region (boolean mask).

    The non-vascular component (:func:`vascular_signal`) overlapping a
    ``FAZ_SEED_MM`` disk at the image center, the largest if several, holes
    filled. Its area is the OCTA literature's FAZ area and is kept as a
    diagnostic; the scored statistic is the leak-proof clear radius
    (:func:`faz_metrics`). Both read the angiogram, not a label mask:
    expert OCTA labels omit perifoveal capillaries, so a zone read from
    them leaks into the intercapillary spaces (0.9 mm2 on ROSE-1 against
    0.37 from the angiograms, and no closing radius repairs it).
    """
    labels, _ = ndimage.label(~vascular_signal(signal, mm_per_px))
    rows, cols = np.indices(signal.shape)
    from_center = np.hypot(rows - signal.shape[0] / 2, cols - signal.shape[1] / 2)
    seed = (from_center <= FAZ_SEED_MM / mm_per_px) & (labels > 0)
    if not seed.any():
        return np.zeros(signal.shape, dtype=bool)
    candidates = np.unique(labels[seed])
    sizes = ndimage.sum(np.ones_like(labels), labels, index=candidates)
    return ndimage.binary_fill_holes(labels == candidates[np.argmax(sizes)])


def faz_metrics(signal: np.ndarray, mm_per_px: float) -> dict[str, float]:
    """FAZ clear radius (mm, scored), with the zone's area (mm2) and acircularity as diagnostics."""
    radius_px = clear_radius_px(vascular_signal(signal, mm_per_px), FAZ_SEARCH_MM / mm_per_px)
    zone = avascular_zone(signal, mm_per_px)
    area_px = int(zone.sum())
    return {
        "faz_radius_mm": radius_px * mm_per_px,
        "faz_area_mm2": area_px * mm_per_px**2,
        "faz_acircularity": (
            float(perimeter(zone, neighborhood=8) ** 2 / (4 * np.pi * area_px))
            if area_px
            else float("nan")
        ),
    }


def octa_window(
    edges: pd.DataFrame, fovea_center: tuple[float, float], radii: np.ndarray | None = None
) -> np.ndarray:
    """The network as an OCTA en-face scan: a 3 x 3 mm window on the fovea at ROSE scale.

    OCTA images flow, not caliber: an 8 um capillary is as bright as a
    pixel, so every segment is drawn at least one pixel wide here, unlike
    the fundus raster, where vessels below half a pixel vanish.
    """
    units_per_px = OCTA_SCAN_MM / MM_PER_UNIT / OCTA_SIZE_PX
    if radii is not None:
        radii = np.maximum(np.asarray(radii, dtype=float), 0.6 * units_per_px)
    return rasterize_window(
        edges, fovea_center, OCTA_SCAN_MM / MM_PER_UNIT, OCTA_SIZE_PX, radii
    )


def capillary_statistics(
    skeleton: np.ndarray, excluded: np.ndarray, mm_per_px: float
) -> dict[str, float]:
    """Capillary-scale spacing and length density of a vessel skeleton.

    Over the pixels outside ``excluded`` (the FAZ, whose emptiness is a
    different statistic): ``octa_intervessel_um`` is twice the mean distance
    from a non-vessel pixel to the nearest skeleton pixel (parallel vessels
    a distance s apart read s/2), and ``octa_skeleton_mm_per_mm2`` is
    skeleton length per tissue area. On ROSE-1 the skeleton is the expert
    label (large vessels and capillary centerlines); on the model it is the
    skeleton of :func:`octa_window`. The labels omit some perifoveal
    capillaries, so the reference spacing reads a little wide and the
    density a little low; the model is measured without that omission.
    """
    skeleton = np.asarray(skeleton, dtype=bool)
    outside = ~np.asarray(excluded, dtype=bool)
    if not skeleton.any() or not outside.any():
        return {"octa_intervessel_um": float("nan"), "octa_skeleton_mm_per_mm2": 0.0}
    clearance = ndimage.distance_transform_edt(~skeleton)
    gaps = outside & ~skeleton
    return {
        "octa_intervessel_um": float(2 * clearance[gaps].mean() * mm_per_px * 1000),
        "octa_skeleton_mm_per_mm2": float(
            skeleton[outside].sum() * mm_per_px / (outside.sum() * mm_per_px**2)
        ),
    }


def site_reach(
    pop: pd.DataFrame, sites: np.ndarray, perfusion_radius: float, vessel_type: int | None
) -> np.ndarray:
    """Mask of demand sites within perfusion_radius of a frozen (typed) vessel."""
    frozen = pop[pop.frozen]
    if vessel_type is not None:
        frozen = frozen[frozen.vessel_type == vessel_type]
    positions = frozen[["x", "y", "z"]].to_numpy(dtype=float)
    if len(positions) == 0 or len(sites) == 0:
        return np.zeros(len(sites), dtype=bool)
    distances, _ = cKDTree(positions).query(sites, k=1)
    return distances <= perfusion_radius


def perfused_fraction(
    pop: pd.DataFrame,
    semi_axes,
    site_spacing: float,
    perfusion_radius: float,
    vessel_type: int | None = None,
) -> float:
    """Fraction of tissue demand sites within perfusion_radius of a frozen vessel.

    The direct measure of how completely the network colonizes its territory
    (roadmap idea 2). Mirrors the PerfusionDemand component's site lattice.
    With ``vessel_type`` given, only frozen vessels of that type count (e.g.
    arterial supply coverage vs. venous drainage coverage). This is the
    growth-completeness measure; :func:`paired_perfused_fraction` is the
    physiological one.
    """
    from vivarium_eye_vessels.components.boundaries import generate_demand_sites

    sites = generate_demand_sites(np.asarray(semi_axes, dtype=float), site_spacing)
    if len(sites) == 0:
        return 0.0
    return float(site_reach(pop, sites, perfusion_radius, vessel_type).mean())


def paired_perfused_fraction(
    pop: pd.DataFrame, semi_axes, site_spacing: float, perfusion_radius: float
) -> float:
    """Fraction of demand sites within perfusion_radius of both an artery and a vein.

    Tissue is perfused when blood can arrive and leave: a capillary bed needs
    an arteriole to feed it and a venule to drain it. A site reached by one
    tree only is colonized, not perfused, so this is the scored perfusion
    target; the any-vessel fraction measures growth completeness.
    """
    from vivarium_eye_vessels.components.boundaries import generate_demand_sites

    sites = generate_demand_sites(np.asarray(semi_axes, dtype=float), site_spacing)
    if len(sites) == 0:
        return 0.0
    supplied = site_reach(pop, sites, perfusion_radius, VESSEL_TYPE_ARTERY)
    drained = site_reach(pop, sites, perfusion_radius, VESSEL_TYPE_VEIN)
    return float((supplied & drained).mean())


def arcade_caliber_ratio(
    pop: pd.DataFrame, disc_center, zone: tuple[float, float] = AVR_ZONE
) -> float:
    """Mean artery over mean vein caliber of the depth-0 arcades near the disc.

    Clinical AVR (CRAE/CRVE) is read on the major arcades in a fixed zone
    around the disc margin. Averaging each trunk's whole length instead
    makes the ratio depend on how far each tapering trunk happens to run
    -- a seed lottery of 0.63-0.85 at a configured 0.67 -- so only depth-0
    particles between ``zone[0]`` and ``zone[1]`` from the disc center
    count, and each trunk (depth-0 path) counts once, as each vessel does
    in CRAE/CRVE: a trunk that coils inside the zone while tapering to a
    thread would otherwise swamp the particle mean (1.40 on one sweep
    seed, 210 points of a 242-point score). NaN when either tree is
    absent from the zone.
    """
    trunks = pop[(pop.depth == 0) & (pop.radius > 0) & (pop.vessel_type > 0)]
    distance = np.hypot(trunks.x - disc_center[0], trunks.y - disc_center[1])
    in_zone = trunks[(distance >= zone[0]) & (distance < zone[1])]
    per_trunk = in_zone.groupby(["vessel_type", "path_id"]).radius.mean()
    by_type = per_trunk.groupby(level="vessel_type").mean()
    if VESSEL_TYPE_ARTERY not in by_type.index or VESSEL_TYPE_VEIN not in by_type.index:
        return float("nan")
    return float(by_type[VESSEL_TYPE_ARTERY] / by_type[VESSEL_TYPE_VEIN])


def binarize_mask(mask: np.ndarray, size: int = RASTER_SIZE) -> np.ndarray:
    """Downsample a real vessel mask to the common raster resolution.

    Block averaging followed by a majority vote: a downsampled pixel is
    vessel when more than half of its block is. This keeps calibers true
    (a lower threshold would thicken every thin vessel by a pixel or two)
    and matches the rule :func:`rasterize_network` draws the simulation by.
    """
    mask = np.asarray(mask, dtype=float)
    if mask.max() > 1:
        mask = mask / 255.0
    factor = max(int(np.ceil(max(mask.shape) / size)), 1)
    if factor > 1:
        trim_r = mask.shape[0] - mask.shape[0] % factor
        trim_c = mask.shape[1] - mask.shape[1] % factor
        mask = mask[:trim_r, :trim_c]
        mask = mask.reshape(
            mask.shape[0] // factor, factor, mask.shape[1] // factor, factor
        ).mean(axis=(1, 3))
    return mask > 0.5


def imaged_region(binary: np.ndarray) -> np.ndarray:
    """The convex hull of the vessel pixels: the part of the frame that was imaged."""
    if binary.sum() < 3:  # no hull to speak of
        return binary.astype(bool)
    return convex_hull_image(binary)


########################
# Image-based metrics  #
########################


def box_counting_dimension(image: np.ndarray, sizes: np.ndarray = BOX_SIZES) -> float:
    """Box-counting fractal dimension of a binary image over a fixed box range."""
    if image.sum() < 10:
        return float("nan")

    counts = []
    for box in sizes:
        binned = np.add.reduceat(
            np.add.reduceat(image, np.arange(0, image.shape[0], box), axis=0),
            np.arange(0, image.shape[1], box),
            axis=1,
        )
        counts.append((binned > 0).sum())

    coeffs = np.polyfit(np.log(sizes), np.log(counts), 1)
    return float(-coeffs[0])


def neighbor_counts(skeleton: np.ndarray) -> np.ndarray:
    """Number of 8-connected skeleton neighbors of every pixel."""
    return ndimage.convolve(skeleton.astype(int), NEIGHBOR_KERNEL, mode="constant", cval=0)


def prune_spurs(skeleton: np.ndarray, max_length: int = MIN_BRANCH_PIXELS) -> np.ndarray:
    """Remove terminal branches shorter than ``max_length`` pixels.

    Thinning a vessel with a ragged edge leaves short spurs at the edge
    bumps; each one adds a junction that is not a branch point. A spur is a
    branch that ends in an endpoint, touches a junction, and is shorter
    than a countable branch. Removal can expose a new spur, so the pass
    repeats a few times; the skeleton is re-thinned after each pass to
    tidy the junction remnants.
    """
    skeleton = skeleton.astype(bool).copy()
    for _ in range(3):
        counts = neighbor_counts(skeleton)
        junctions = skeleton & (counts >= 3)
        endpoints = skeleton & (counts == 1)
        labels, n_labels = ndimage.label(skeleton & ~junctions, structure=np.ones((3, 3)))
        if n_labels == 0:
            break
        index = np.arange(1, n_labels + 1)
        sizes = ndimage.sum(np.ones_like(labels), labels, index=index)
        has_endpoint = ndimage.maximum(endpoints, labels, index=index) > 0
        beside_junction = neighbor_counts(junctions) > 0
        touches_junction = ndimage.maximum(beside_junction, labels, index=index) > 0
        spur = (sizes < max_length) & has_endpoint & touches_junction
        if not spur.any():
            break
        skeleton &= ~np.isin(labels, index[spur])
        skeleton = skeletonize(skeleton)
    return skeleton


def vessel_skeleton(binary: np.ndarray) -> np.ndarray:
    """The measured centerline: thinned, then pruned of spurs."""
    return prune_spurs(skeletonize(binary))


def chain_arc_length(coords: np.ndarray) -> float:
    """Arc length of a set of 8-connected pixels, in px.

    The minimum spanning tree over neighboring pixels (orthogonal steps 1,
    diagonal steps sqrt 2) — a pixel count would read a 45-degree line as
    30% shorter than it is.
    """
    n_pixels = len(coords)
    if n_pixels < 2:
        return 0.0
    pairs = cKDTree(coords).query_pairs(r=1.5, output_type="ndarray")
    if len(pairs) == 0:
        return 0.0
    weights = np.linalg.norm(coords[pairs[:, 0]] - coords[pairs[:, 1]], axis=1)
    graph = coo_matrix((weights, (pairs[:, 0], pairs[:, 1])), shape=(n_pixels, n_pixels))
    return float(minimum_spanning_tree(graph).sum())


def skeleton_branches(skeleton: np.ndarray, binary: np.ndarray | None = None) -> pd.DataFrame:
    """Decompose a skeleton into branches between junctions.

    Junction pixels (3+ skeleton neighbors) are removed; each remaining
    connected component of at least ``MIN_BRANCH_PIXELS`` is one branch.
    Returns per-branch arc lengths and endpoint chord distances (their ratio
    is the tortuosity). When the ``binary`` vessel image is given, each
    branch also gets a ``diameter_px``: twice the mean distance-transform
    value along the branch (the medial-axis width estimate), which recovers
    local caliber even from masks that carry no explicit radii.
    """
    skeleton = skeleton.astype(bool)
    junctions = skeleton & (neighbor_counts(skeleton) >= 3)
    branches = skeleton & ~junctions

    edt = None if binary is None else ndimage.distance_transform_edt(binary)
    labels, n_labels = ndimage.label(branches, structure=np.ones((3, 3)))
    records = []
    for slc, label in zip(ndimage.find_objects(labels), range(1, n_labels + 1)):
        branch_img = labels[slc] == label
        coords = np.argwhere(branch_img)
        if len(coords) < MIN_BRANCH_PIXELS:
            continue
        # Endpoints: branch pixels with at most one neighbor within the branch
        ends = np.argwhere(branch_img & (neighbor_counts(branch_img) <= 1))
        if len(ends) >= 2:
            dists = np.linalg.norm(ends[:, None, :] - ends[None, :, :], axis=-1)
            chord = dists.max()
        else:
            chord = float(np.linalg.norm(coords.max(axis=0) - coords.min(axis=0)))
        if chord < 1:
            continue
        record = {"length_px": chain_arc_length(coords.astype(float)), "chord_px": chord}
        if edt is not None:
            offset = np.array([s.start for s in slc])
            rows, cols = (coords + offset).T
            record["diameter_px"] = float(2.0 * edt[rows, cols].mean())
        records.append(record)

    columns = ["length_px", "chord_px"] + ([] if edt is None else ["diameter_px"])
    frame = pd.DataFrame(records, columns=columns)
    frame["tortuosity"] = frame.length_px / frame.chord_px
    return frame


def skeleton_pixel_diameters(
    binary: np.ndarray, skeleton: np.ndarray | None = None
) -> np.ndarray:
    """Local vessel diameter (2 x EDT) at every skeleton pixel, in px.

    Each skeleton pixel is one unit of centerline length at its local
    caliber, so this distribution is the length-weighted caliber profile
    of the network — matching it matches "length x width" without any
    binning into diameter strata.
    """
    if skeleton is None:
        skeleton = vessel_skeleton(binary)
    edt = ndimage.distance_transform_edt(binary)
    return 2.0 * edt[skeleton]


def wide_junction_spacing(
    skeleton: np.ndarray, binary: np.ndarray, min_diameter_px: float = 4.0
) -> float:
    """Mean skeleton distance between junctions along wide vessels, in px.

    Real arcades throw off side branches at short, regular intervals (a
    comb-like pattern), so the wide-vessel skeleton carries many junctions
    per unit length. Wide pixels are those whose medial-axis diameter
    (2 x EDT) exceeds ``min_diameter_px``; the spacing is wide skeleton
    pixels per junction on a wide pixel. Crossing vessels create the same
    false junctions in simulated rasters and real masks, so the estimator
    stays comparable. NaN when there is no wide skeleton.
    """
    skeleton = skeleton.astype(bool)
    junctions = skeleton & (neighbor_counts(skeleton) >= 3)
    edt = ndimage.distance_transform_edt(binary)
    wide = skeleton & (2.0 * edt > min_diameter_px)
    n_wide = int(wide.sum())
    if n_wide == 0:
        return float("nan")
    # One skeleton joint spans several junction pixels; count connected
    # junction clusters so a branch point counts once
    labels, n_clusters = ndimage.label(junctions, structure=np.ones((3, 3)))
    on_wide = np.unique(labels[wide & junctions])
    n_junctions = int((on_wide > 0).sum())
    return float(n_wide / max(n_junctions, 1))


####################
# Arcade geometry  #
####################

WIDE_DIAMETER_PX = 4.0  # the arcade class: wide_share, junction spacing, reach
THICK_DIAMETER_PX = 6.0  # the caliber a real arcade rarely exceeds at this raster
ARCADE_MIN_DISC_DISTANCE_PX = 100.0  # alignment is read clear of the disc's convergence


def local_orientation(skeleton: np.ndarray, size: int = 7) -> np.ndarray:
    """Unit tangent (d_row, d_col) of the skeleton at every pixel.

    The principal axis of the skeleton pixels in a ``size`` x ``size`` window,
    so a vessel's direction is read over a few pixels of centerline rather
    than one pixel step. Sign is arbitrary (a tangent line, not a heading).
    """
    mask = skeleton.astype(float)
    rows, cols = np.indices(skeleton.shape, dtype=float)
    count = np.maximum(ndimage.uniform_filter(mask, size, mode="constant"), 1e-12)

    def local_mean(field: np.ndarray) -> np.ndarray:
        return ndimage.uniform_filter(mask * field, size, mode="constant") / count

    mean_r, mean_c = local_mean(rows), local_mean(cols)
    var_r = local_mean(rows * rows) - mean_r**2
    var_c = local_mean(cols * cols) - mean_c**2
    cov = local_mean(rows * cols) - mean_r * mean_c
    theta = 0.5 * np.arctan2(2 * cov, var_r - var_c)
    return np.stack([np.cos(theta), np.sin(theta)], axis=-1)


def disc_center(
    coords: np.ndarray, tangents: np.ndarray, weights: np.ndarray, iterations: int = 5
) -> np.ndarray:
    """Where the arcades converge: the point nearest the wide vessels' tangent lines.

    Every major vessel radiates from the optic disc, so the disc is the
    point that minimizes the (weighted) squared distance to the tangent
    lines through the wide skeleton pixels — the same estimate for a real
    mask, whose disc is unknown, and a simulated raster, whose disc is
    known (it lands within ~40 px of it). Side branches leave the arcades
    at right angles and their lines miss the disc, so the least squares is
    re-weighted a few times to discount lines that pass far from the
    current estimate (Cauchy weights, 60 px scale).
    """
    normals = np.column_stack([-tangents[:, 1], tangents[:, 0]])
    offsets = (normals * coords).sum(axis=1)
    current = weights.astype(float)
    center = coords.mean(axis=0)
    for _ in range(iterations):
        system = np.einsum("i,ij,ik->jk", current, normals, normals)
        rhs = np.einsum("i,ij,i->j", current, normals, offsets)
        center = np.linalg.solve(system, rhs)
        miss = np.abs(normals @ center - offsets)
        current = weights / (1.0 + (miss / 60.0) ** 2)
    return center


def arcade_geometry(binary: np.ndarray, disc: np.ndarray | None = None) -> dict[str, float]:
    """Long-scale geometry of the arcades in a fundus-sized binary image.

    Real arcades leave the disc and run to the periphery, straight at the
    scale of the image and tapering as they go; branch-level tortuosity
    (a few pixels between junctions) cannot see a trunk that curls back
    on itself over 100 px. Three statistics read that scale, all on the
    wide (> ``WIDE_DIAMETER_PX``) skeleton and relative to the disc, which
    is estimated from the image itself (:func:`disc_center`) unless given:

    - ``arcade_radial_alignment``: mean |cos| between the local tangent and
      the direction from the disc, over wide pixels farther than
      ``ARCADE_MIN_DISC_DISTANCE_PX`` from it (1 = every arcade points
      away from the disc, 0 = rings around it)
    - ``arcade_reach_px``: mean distance of wide skeleton pixels from the
      disc — how far the arcade class extends into the field
    - ``thick_share``: share of skeleton length wider than
      ``THICK_DIAMETER_PX`` — the top of the caliber profile, where the KS
      statistic has too little mass to notice a trunk twice too wide

    NaN entries when the image has no wide skeleton.
    """
    skeleton = vessel_skeleton(binary)
    diameters = 2.0 * ndimage.distance_transform_edt(binary)
    wide = skeleton & (diameters > WIDE_DIAMETER_PX)
    if wide.sum() < 3:
        return {
            "disc_row_px": float("nan"),
            "disc_col_px": float("nan"),
            "arcade_radial_alignment": float("nan"),
            "arcade_reach_px": float("nan"),
            "thick_share": float((diameters[skeleton] > THICK_DIAMETER_PX).mean()),
        }
    coords = np.argwhere(wide).astype(float)
    tangents = local_orientation(skeleton)[wide]
    if disc is None:
        disc = disc_center(coords, tangents, diameters[wide])
    radial = coords - disc
    distance = np.linalg.norm(radial, axis=1)
    far = distance > ARCADE_MIN_DISC_DISTANCE_PX
    alignment = np.abs((tangents[far] * radial[far]).sum(axis=1) / distance[far])
    return {
        "disc_row_px": float(disc[0]),
        "disc_col_px": float(disc[1]),
        "arcade_radial_alignment": float(alignment.mean()) if far.any() else float("nan"),
        "arcade_reach_px": float(distance.mean()),
        "thick_share": float((diameters[skeleton] > THICK_DIAMETER_PX).mean()),
    }


def image_metrics(binary: np.ndarray) -> dict[str, Any]:
    """All image-based metrics for one binary vessel image.

    Densities are fractions of the imaged region (the vessel pixels' convex
    hull), whose share of the frame is reported as ``fov_fraction``. The
    arcade geometry (:func:`arcade_geometry`) is read over the whole image,
    which for a real mask is the fundus field; a simulated raster larger
    than that is windowed first (:func:`fundus_window`).
    """
    skeleton = vessel_skeleton(binary)
    region = imaged_region(binary)
    branches = skeleton_branches(skeleton, binary)
    return {
        "fov_fraction": float(region.mean()),
        "fractal_dimension": box_counting_dimension(skeleton),
        "skeleton_density": float(skeleton[region].mean()) if region.any() else 0.0,
        "area_density": float(binary[region].mean()) if region.any() else 0.0,
        "n_branches": int(len(branches)),
        "branch_length_px": branches.length_px.tolist(),
        "branch_tortuosity": branches.tortuosity.tolist(),
        "branch_diameter_px": branches.diameter_px.tolist(),
        "pixel_diameter_px": skeleton_pixel_diameters(binary, skeleton).tolist(),
        "wide_junction_spacing_px": wide_junction_spacing(skeleton, binary),
        "macular_clear_radius_px": clear_radius_px(
            skeleton, MACULA_SEARCH_MM / FUNDUS_MM_PER_PX
        ),
        **arcade_geometry(binary),
    }


#######################
# Tree-based metrics  #
#######################


CONTINUATION_CALIBER_RATIO = 0.98


def true_bifurcations(pop: pd.DataFrame):
    """Yield (parent, daughters) at true bifurcations of the vessel tree.

    A parent's continuation child (created by PathFreezer, at nearly the
    parent's radius) is not a daughter branch, so it is excluded; only
    parents with two or more daughters count. The continuation is recognized
    by caliber rather than by path label: it is the one child that keeps at
    least ``CONTINUATION_CALIBER_RATIO`` of the parent's radius, while split
    daughters are born at Murray calibers below that. A freezer continuation
    never comes in pairs, so when two or more children keep the caliber
    (daughters clipped to the caliber floor, or an uncalibered tree) they are
    all daughters. Sprout points — a continuation plus one side branch — are
    therefore skipped rather than fit as if they were bifurcations.
    """
    children = pop[pop.parent_id >= 0]
    children = children[children.parent_id.isin(pop.index)]
    calibered = "radius" in pop.columns
    for parent_id, group in children.groupby("parent_id"):
        parent = pop.loc[parent_id]
        daughters = group
        if calibered:
            keeps_caliber = group.radius >= CONTINUATION_CALIBER_RATIO * parent.radius
            if keeps_caliber.sum() == 1:
                daughters = group[~keeps_caliber]
        if len(daughters) >= 2:
            yield parent, daughters


def bifurcation_angles(pop: pd.DataFrame) -> np.ndarray:
    """Angles (degrees) between daughter segments at true tree bifurcations."""
    angles = []
    for parent, group in true_bifurcations(pop):
        vectors = group[["x", "y", "z"]].values - parent[["x", "y", "z"]].values.astype(float)
        norms = np.linalg.norm(vectors, axis=1)
        valid = norms > 1e-12
        vectors, norms = vectors[valid], norms[valid]
        if len(vectors) < 2:
            continue
        v1, v2 = vectors[0] / norms[0], vectors[1] / norms[1]
        cos_angle = np.clip(np.dot(v1, v2), -1.0, 1.0)
        angles.append(np.degrees(np.arccos(cos_angle)))
    return np.array(angles)


def junction_exponents(pop: pd.DataFrame) -> np.ndarray:
    """Fitted junction exponents k with r0**k == r1**k + r2**k at true bifurcations.

    Murray's law predicts k close to 3. Bifurcations without calibers, or
    where a daughter is at least as wide as the parent (e.g. from caliber
    flooring), are skipped since no exponent exists there.
    """
    exponents = []
    for parent, group in true_bifurcations(pop):
        r0 = float(parent.radius)
        r1, r2 = (float(r) for r in group.radius.values[:2])
        if min(r0, r1, r2) <= 0 or max(r1, r2) >= r0:
            continue

        def mismatch(k: float) -> float:
            return (r1 / r0) ** k + (r2 / r0) ** k - 1.0

        if mismatch(0.5) <= 0 or mismatch(15.0) >= 0:
            continue
        exponents.append(brentq(mismatch, 0.5, 15.0))
    return np.array(exponents)


def graph_cycles(pop: pd.DataFrame) -> int:
    """Number of independent cycles in the vessel graph (E - N + C).

    Parent-child edges alone form a forest (zero cycles); every anastomosis
    that joins two already-connected components adds exactly one cycle, so
    this counts the loops that make the network perfusable rather than
    tree-like.
    """
    in_table = pop.parent_id.isin(pop.index)
    edges = list(zip(pop.index[in_table], pop.parent_id[in_table]))
    if "anastomosis_id" in pop.columns:
        joined = pop.anastomosis_id.isin(pop.index) & (pop.anastomosis_id >= 0)
        edges.extend(zip(pop.index[joined], pop.anastomosis_id[joined]))
    if not edges:
        return 0

    nodes = {node for edge in edges for node in edge}
    parent = {node: node for node in nodes}

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    components = len(nodes)
    for a, b in edges:
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[root_a] = root_b
            components -= 1
    return len(edges) - len(nodes) + components


def tree_segment_lengths(pop: pd.DataFrame) -> np.ndarray:
    """Arc lengths of vessel segments between consecutive branch points.

    A segment's lower end is a leaf or a branch point (a node with two or
    more children of any kind, side branches included); its length
    accumulates every tree edge walking up until the next branch point or
    the root. Unlike the raster-based branch lengths, these are measured in
    simulation units directly on the tree, independent of projection,
    caliber, and raster resolution.
    """
    on_path = pop[pop.path_id >= 0]
    members = on_path.index
    parents = on_path.parent_id.where(on_path.parent_id.isin(members), -1)
    child_counts = parents[parents >= 0].value_counts()
    parent_map = parents.to_dict()
    coords = {idx: row for idx, row in zip(members, on_path[["x", "y", "z"]].to_numpy(float))}

    lengths = []
    for start in members:
        if child_counts.get(start, 0) == 1:
            continue  # interior continuation node, not a segment's lower end
        length = 0.0
        node = start
        while True:
            parent = parent_map.get(node, -1)
            if parent < 0:
                break
            length += float(np.linalg.norm(coords[node] - coords[parent]))
            node = parent
            if child_counts.get(node, 0) >= 2:
                break
        if length > 0:
            lengths.append(length)
    return np.array(lengths)


def path_chains(pop: pd.DataFrame):
    """Yield the root-to-leaf coordinate array of every chain of a path.

    A ``path_id`` group is a chain of particles linked by ``parent_id``
    (occasionally two chains, since sibling branches from one split share a
    path id); each leaf is walked back to the group's root.
    """
    on_path = pop[pop.path_id >= 0]
    for _, group in on_path.groupby("path_id"):
        members = set(group.index)
        parents_in_group = set(group.parent_id) & members
        leaves = [idx for idx in group.index if idx not in parents_in_group]
        for leaf in leaves:
            chain = [leaf]
            while True:
                parent = group.parent_id.get(chain[-1], -1)
                if parent in members:
                    chain.append(parent)
                else:
                    break
            yield group.loc[chain[::-1], ["x", "y", "z"]].to_numpy(dtype=float)


def path_tortuosity(pop: pd.DataFrame) -> np.ndarray:
    """Arc-length over chord-length for each root-to-leaf chain of a path."""
    ratios = []
    for points in path_chains(pop):
        if len(points) < 3:
            continue
        arc = np.linalg.norm(np.diff(points, axis=0), axis=1).sum()
        chord = np.linalg.norm(points[-1] - points[0])
        if chord > 1e-12:
            ratios.append(arc / chord)
    return np.array(ratios)


def path_turning_coherence(pop: pd.DataFrame, min_turns: int = 6) -> np.ndarray:
    """Lag-1 autocorrelation of the signed in-plane turning angle along each chain.

    Consecutive turns of the same sign (an arc) score toward +1; turns that
    flip sign at every node (jitter) score toward -1. This separates
    coherent curvature from step-scale wiggle, which arc-over-chord alone
    cannot, so it is the direct readout of the steering persistence
    (roadmap idea 7). Chains with fewer than ``min_turns`` turns are skipped.
    """
    coherence = []
    for points in path_chains(pop):
        steps = np.diff(points[:, :2], axis=0)
        headings = np.arctan2(steps[:, 1], steps[:, 0])
        turns = np.angle(np.exp(1j * np.diff(headings)))
        if len(turns) < min_turns or np.std(turns) < 1e-9:
            continue  # too short, or constant curvature (autocorrelation undefined)
        coherence.append(np.corrcoef(turns[:-1], turns[1:])[0, 1])
    return np.array(coherence)


#############
# Summaries #
#############


def summarize(values) -> dict[str, float]:
    """Median and quartiles of a distribution, NaN-safe."""
    values = np.asarray(list(values), dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"n": 0, "median": float("nan"), "q25": float("nan"), "q75": float("nan")}
    return {
        "n": int(len(values)),
        "median": float(np.median(values)),
        "q25": float(np.quantile(values, 0.25)),
        "q75": float(np.quantile(values, 0.75)),
    }

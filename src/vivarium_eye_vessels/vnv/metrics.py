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
of the field of view. They support before/after comparison across model
versions, not absolute physiological claims.
"""

from typing import Any

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.optimize import brentq
from scipy.spatial import cKDTree
from skimage.draw import line as draw_line
from skimage.morphology import dilation, disk, skeletonize

RASTER_SIZE = 1024
MIN_BRANCH_PIXELS = 5


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
        ``edges``). When given, segments are drawn with their caliber instead
        of one pixel wide.
    """
    a, b = bounds
    aspect = b / a
    width = size if a >= b else max(int(size * a / b), 1)
    height = max(int(width * aspect), 1)

    image = np.zeros((height, width), dtype=bool)
    if edges.empty:
        return image

    def to_pixel(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        col = np.clip(((x + a) / (2 * a) * (width - 1)).round().astype(int), 0, width - 1)
        row = np.clip(((y + b) / (2 * b) * (height - 1)).round().astype(int), 0, height - 1)
        return row, col

    r0, c0 = to_pixel(edges.x0.values, edges.y0.values)
    r1, c1 = to_pixel(edges.x1.values, edges.y1.values)

    if radii is None:
        widths_px = np.ones(len(edges), dtype=int)
    else:
        px_per_unit = (width - 1) / (2 * a)
        widths_px = np.maximum(
            np.round(2 * np.asarray(radii, dtype=float) * px_per_unit).astype(int), 1
        )

    # Draw segments in buckets of equal pixel width, thickened by dilation
    for width_px in np.unique(widths_px):
        layer = np.zeros_like(image)
        for i in np.nonzero(widths_px == width_px)[0]:
            rr, cc = draw_line(r0[i], c0[i], r1[i], c1[i])
            layer[rr, cc] = True
        dilation_radius = int(round((width_px - 1) / 2))
        if dilation_radius > 0:
            layer = dilation(layer, disk(dilation_radius))
        image |= layer
    return image


def perfused_fraction(
    pop: pd.DataFrame,
    semi_axes,
    site_spacing: float,
    perfusion_radius: float,
) -> float:
    """Fraction of tissue demand sites within perfusion_radius of a frozen vessel.

    The direct measure of how completely the network colonizes its territory
    (roadmap idea 2). Mirrors the PerfusionDemand component's site lattice.
    """
    from vivarium_eye_vessels.components.boundaries import generate_demand_sites

    sites = generate_demand_sites(np.asarray(semi_axes, dtype=float), site_spacing)
    frozen = pop[pop.frozen][["x", "y", "z"]].to_numpy(dtype=float)
    if len(frozen) == 0 or len(sites) == 0:
        return 0.0
    distances, _ = cKDTree(frozen).query(sites, k=1)
    return float((distances <= perfusion_radius).mean())


def binarize_mask(mask: np.ndarray, size: int = RASTER_SIZE) -> np.ndarray:
    """Downsample a real vessel mask to the common raster resolution.

    Uses block averaging followed by a low threshold so that thin vessels
    survive the downsampling.
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
    return mask > 0.1


########################
# Image-based metrics  #
########################


def box_counting_dimension(image: np.ndarray) -> float:
    """Box-counting fractal dimension of a binary image."""
    pixels = np.argwhere(image)
    if len(pixels) < 10:
        return float("nan")

    max_dim = max(image.shape)
    max_exp = int(np.floor(np.log2(max_dim / 4)))
    sizes = 2 ** np.arange(1, max_exp + 1)

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


def skeleton_branches(skeleton: np.ndarray) -> pd.DataFrame:
    """Decompose a skeleton into branches between junctions.

    Junction pixels (3+ skeleton neighbors) are removed; each remaining
    connected component is one branch. Returns per-branch pixel counts
    (length) and endpoint chord distances (for tortuosity).
    """
    skeleton = skeleton.astype(bool)
    kernel = np.ones((3, 3), dtype=int)
    kernel[1, 1] = 0
    neighbor_count = ndimage.convolve(skeleton.astype(int), kernel, mode="constant", cval=0)
    junctions = skeleton & (neighbor_count >= 3)
    branches = skeleton & ~junctions

    labels, n_labels = ndimage.label(branches, structure=np.ones((3, 3)))
    records = []
    for slc, label in zip(ndimage.find_objects(labels), range(1, n_labels + 1)):
        coords = np.argwhere(labels[slc] == label)
        n_pixels = len(coords)
        if n_pixels < MIN_BRANCH_PIXELS:
            continue
        # Endpoints: branch pixels with at most one neighbor within the branch
        branch_img = labels[slc] == label
        local_neighbors = ndimage.convolve(
            branch_img.astype(int), kernel, mode="constant", cval=0
        )
        ends = np.argwhere(branch_img & (local_neighbors <= 1))
        if len(ends) >= 2:
            dists = np.linalg.norm(ends[:, None, :] - ends[None, :, :], axis=-1)
            chord = dists.max()
        else:
            chord = float(np.linalg.norm(coords.max(axis=0) - coords.min(axis=0)))
        if chord < 1:
            continue
        records.append({"length_px": n_pixels, "chord_px": chord})

    frame = pd.DataFrame(records, columns=["length_px", "chord_px"])
    frame["tortuosity"] = frame.length_px / frame.chord_px
    return frame


def image_metrics(binary: np.ndarray) -> dict[str, Any]:
    """All image-based metrics for one binary vessel image."""
    skeleton = skeletonize(binary)
    branches = skeleton_branches(skeleton)
    return {
        "fractal_dimension": box_counting_dimension(skeleton),
        "skeleton_density": float(skeleton.mean()),
        "area_density": float(binary.mean()),
        "n_branches": int(len(branches)),
        "branch_length_px": branches.length_px.tolist(),
        "branch_tortuosity": branches.tortuosity.tolist(),
    }


#######################
# Tree-based metrics  #
#######################


def bifurcation_angles(pop: pd.DataFrame) -> np.ndarray:
    """Angles (degrees) between daughter segments at tree bifurcations."""
    children = pop[pop.parent_id >= 0]
    children = children[children.parent_id.isin(pop.index)]
    angles = []
    for parent_id, group in children.groupby("parent_id"):
        if len(group) < 2:
            continue
        parent = pop.loc[parent_id]
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
    """Fitted junction exponents k with r0**k == r1**k + r2**k at bifurcations.

    Murray's law predicts k close to 3. Bifurcations without calibers, or
    where a daughter is at least as wide as the parent (e.g. from caliber
    flooring), are skipped since no exponent exists there.
    """
    children = pop[pop.parent_id >= 0]
    children = children[children.parent_id.isin(pop.index)]
    exponents = []
    for parent_id, group in children.groupby("parent_id"):
        if len(group) < 2:
            continue
        r0 = float(pop.loc[parent_id, "radius"])
        r1, r2 = (float(r) for r in group.radius.values[:2])
        if min(r0, r1, r2) <= 0 or max(r1, r2) >= r0:
            continue

        def mismatch(k: float) -> float:
            return (r1 / r0) ** k + (r2 / r0) ** k - 1.0

        if mismatch(0.5) <= 0 or mismatch(15.0) >= 0:
            continue
        exponents.append(brentq(mismatch, 0.5, 15.0))
    return np.array(exponents)


def path_tortuosity(pop: pd.DataFrame) -> np.ndarray:
    """Arc-length over chord-length for each root-to-leaf chain of a path.

    A ``path_id`` group is a chain of particles linked by ``parent_id``
    (occasionally two chains, since sibling branches from one split share a
    path id); each leaf is walked back to the group's root.
    """
    ratios = []
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
            if len(chain) < 3:
                continue
            points = group.loc[chain, ["x", "y", "z"]].values
            steps = np.linalg.norm(np.diff(points, axis=0), axis=1)
            arc = steps.sum()
            chord = np.linalg.norm(points[-1] - points[0])
            if chord > 1e-12:
                ratios.append(arc / chord)
    return np.array(ratios)


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

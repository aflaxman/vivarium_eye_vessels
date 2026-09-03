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
from skimage.morphology import convex_hull_image, skeletonize

RASTER_SIZE = 1024
MIN_BRANCH_PIXELS = 5
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
    arterial supply coverage vs. venous drainage coverage).
    """
    from vivarium_eye_vessels.components.boundaries import generate_demand_sites

    sites = generate_demand_sites(np.asarray(semi_axes, dtype=float), site_spacing)
    frozen = pop[pop.frozen]
    if vessel_type is not None:
        frozen = frozen[frozen.vessel_type == vessel_type]
    frozen_positions = frozen[["x", "y", "z"]].to_numpy(dtype=float)
    if len(frozen_positions) == 0 or len(sites) == 0:
        return 0.0
    distances, _ = cKDTree(frozen_positions).query(sites, k=1)
    return float((distances <= perfusion_radius).mean())


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


def image_metrics(binary: np.ndarray) -> dict[str, Any]:
    """All image-based metrics for one binary vessel image.

    Densities are fractions of the imaged region (the vessel pixels' convex
    hull), whose share of the frame is reported as ``fov_fraction``.
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

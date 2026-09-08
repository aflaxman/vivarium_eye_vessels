"""Tips do not pass through vessels of their own plane.

In a fundus the branches of one tree never cross one another: two vessels in
the same plane that meet must fuse or stop, and the only crossings are the
large arteries passing over the large veins, which lie at different depths.
The model's tips are steered by soft repulsion and a strong pull toward
unperfused tissue, so a tip that met a vessel head-on passed through it and
the raster showed a crossing -- one junction in five was a crossing, against
one in eight in the HRF masks (twenty-sixth pass).
"""

from typing import List

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from vivarium import Component
from vivarium.framework.engine import Builder


def segments_cross(p0, p1, q0, q1) -> np.ndarray:
    """Proper 2-D intersection of segments p0->p1 and q0->q1, vectorized over rows."""

    def cross(a, b):
        return a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]

    r, s = p1 - p0, q1 - q0
    denominator = cross(r, s)
    offset = q0 - p0
    with np.errstate(divide="ignore", invalid="ignore"):
        t = cross(offset, s) / denominator
        u = cross(offset, r) / denominator
    return (np.abs(denominator) > 1e-15) & (t > 0) & (t < 1) & (u > 0) & (u < 1)


class CollisionGuard(Component):
    """Stop a tip whose move would cross a frozen vessel of its own plane.

    Each step, every non-capillary tip in the guarded layers is checked
    against the frozen segments of the same layer near its move (x-y
    projection). A pair may cross only when the two vessels belong to
    different trees and both are at least ``crossing_min_radius`` wide --
    an artery passing over a vein, which the retina keeps at different
    depths. Otherwise the tip stays where it was and its path ends there,
    as an extinguished tip's does; the pruner and the capillary bed decide
    what becomes of the branch. :class:`Particle3D` consults the guard as
    it moves the tips.
    """

    CONFIGURATION_DEFAULTS = {
        "collision_guard": {
            "enabled": False,
            "layers": [0],  # plexus layers whose vessels lie in one plane
            # Tips and obstacles at or below this caliber (the capillary bed)
            # are exempt: a bed is a mesh, and its sprouts are not drawn
            "capillary_radius": 0.0,
            # Opposite-tree pairs both at least this wide may cross
            "crossing_min_radius": 0.005,
            # A tip at least this many times wider than the vessel in its way
            # passes over it whatever the trees (a trunk over a twig; the twig
            # is the one that dips). 0 disables
            "crossing_wide_ratio": 0.0,
            # Frozen particles within this of a tip's move are candidate obstacles
            "search_radius": 0.05,
            # What a blocked tip does: "reflect" reverses the part of its
            # move (and velocity) normal to the vessel it would have crossed,
            # as off a wall, so it turns back into open tissue; "slide" keeps
            # only the part along the vessel (contact guidance -- the trail
            # then hugs the vessel and the raster merges the two); "stop"
            # leaves it where it was and ends its path as an extinguished
            # tip's (that halved the visible skeleton) (twenty-sixth pass)
            "mode": "reflect",
        }
    }

    @property
    def required_attributes(self) -> List[str]:
        return [
            "x",
            "y",
            "frozen",
            "parent_id",
            "path_id",
            "radius",
            "vessel_type",
            "layer_id",
        ]

    def setup(self, builder: Builder) -> None:
        self.config = builder.configuration.collision_guard
        self.enabled = bool(self.config.enabled)
        self.layers = [int(layer) for layer in self.config.layers]
        self.capillary_radius = float(self.config.capillary_radius)
        self.crossing_min_radius = float(self.config.crossing_min_radius)
        self.crossing_wide_ratio = float(self.config.crossing_wide_ratio)
        self.search_radius = float(self.config.search_radius)
        self.mode = str(self.config.mode)
        self.total_blocked = 0

    def blocked_tips(
        self, before: pd.DataFrame, after: pd.DataFrame, pop: pd.DataFrame
    ) -> pd.Index:
        """Indices of tips whose move from ``before`` to ``after`` crosses a vessel they may not."""
        return self.blocked_with_tangents(before, after, pop)[0]

    def blocked_with_tangents(
        self, before: pd.DataFrame, after: pd.DataFrame, pop: pd.DataFrame
    ) -> tuple[pd.Index, np.ndarray]:
        """Blocked tips and, for each, the unit direction of the first vessel its move crosses.

        ``before`` and ``after`` carry the tips' x and y (same index);
        ``pop`` is the whole population, with the frozen vessels.
        """
        empty = (pd.Index([]), np.zeros((0, 2)))
        tips = pop.loc[before.index]
        guarded = (tips.path_id >= 0) & tips.layer_id.isin(self.layers)
        if self.capillary_radius > 0:
            guarded &= tips.radius > self.capillary_radius
        tips = tips[guarded]
        if tips.empty:
            return empty
        frozen = pop[pop.frozen & (pop.parent_id >= 0) & pop.layer_id.isin(self.layers)]
        if self.capillary_radius > 0:
            frozen = frozen[frozen.radius > self.capillary_radius]
        frozen = frozen[frozen.parent_id.isin(pop.index)]
        # The fresh trail behind every other tip (its frozen parent to the
        # tip's position before this move) is a vessel too, frozen or not
        active = pop[~pop.frozen & (pop.path_id >= 0) & (pop.parent_id >= 0)]
        active = active[active.layer_id.isin(self.layers) & active.parent_id.isin(pop.index)]
        if self.capillary_radius > 0:
            active = active[active.radius > self.capillary_radius]
        frozen = pd.concat([frozen, active])
        if frozen.empty:
            return empty
        parents = pop.loc[frozen.parent_id]
        obstacle_a = frozen[["x", "y"]].to_numpy(dtype=float)
        obstacle_b = parents[["x", "y"]].to_numpy(dtype=float)

        start = before.loc[tips.index, ["x", "y"]].to_numpy(dtype=float)
        end = after.loc[tips.index, ["x", "y"]].to_numpy(dtype=float)
        tree = cKDTree((obstacle_a + obstacle_b) / 2)
        neighbor_lists = tree.query_ball_point((start + end) / 2, self.search_radius)
        counts = np.fromiter((len(n) for n in neighbor_lists), dtype=int, count=len(tips))
        if counts.sum() == 0:
            return empty
        tip_idx = np.repeat(np.arange(len(tips)), counts)
        obstacle_idx = np.concatenate(
            [np.asarray(n, dtype=int) for n in neighbor_lists if len(n)]
        )

        # A tip's own fresh trail (the segments at its parent, and its own) is not an obstacle
        tip_parent = tips.parent_id.to_numpy()[tip_idx]
        tip_id = tips.index.to_numpy()[tip_idx]
        child_id = frozen.index.to_numpy()[obstacle_idx]
        parent_id = frozen.parent_id.to_numpy()[obstacle_idx]
        own = (child_id == tip_parent) | (parent_id == tip_parent) | (child_id == tip_id)
        crossing = segments_cross(
            start[tip_idx], end[tip_idx], obstacle_a[obstacle_idx], obstacle_b[obstacle_idx]
        )
        crossing &= ~own
        # Large vessels of opposite trees pass over one another
        other_tree = (
            tips.vessel_type.to_numpy()[tip_idx]
            != frozen.vessel_type.to_numpy()[obstacle_idx]
        )
        both_wide = (
            tips.radius.to_numpy(dtype=float)[tip_idx] >= self.crossing_min_radius
        ) & (frozen.radius.to_numpy(dtype=float)[obstacle_idx] >= self.crossing_min_radius)
        allowed = other_tree & both_wide
        if self.crossing_wide_ratio > 0:
            allowed |= tips.radius.to_numpy(dtype=float)[tip_idx] >= (
                self.crossing_wide_ratio * frozen.radius.to_numpy(dtype=float)[obstacle_idx]
            )
        blocked = crossing & ~allowed
        if not blocked.any():
            return empty
        tip_idx, obstacle_idx = tip_idx[blocked], obstacle_idx[blocked]
        # The first vessel along the move decides the tangent
        moves = end[tip_idx] - start[tip_idx]
        along = obstacle_b[obstacle_idx] - obstacle_a[obstacle_idx]
        offset = obstacle_a[obstacle_idx] - start[tip_idx]
        denominator = moves[:, 0] * along[:, 1] - moves[:, 1] * along[:, 0]
        with np.errstate(divide="ignore", invalid="ignore"):
            t_along_move = (
                offset[:, 0] * along[:, 1] - offset[:, 1] * along[:, 0]
            ) / denominator
        order = np.lexsort((t_along_move, tip_idx))
        first = order[np.r_[True, np.diff(tip_idx[order]) != 0]]
        tangents = along[first] / np.maximum(
            np.linalg.norm(along[first], axis=1, keepdims=True), 1e-12
        )
        return tips.index[tip_idx[first]], tangents

    def deflect(
        self, before: pd.DataFrame, updates: pd.DataFrame, pop: pd.DataFrame
    ) -> pd.Index:
        """Apply the guard to a move in place; return the blocked tips.

        ``updates`` holds the moved positions and velocities and is edited:
        in slide mode a blocked tip keeps only the part of its move (and of
        its velocity) along the vessel it would have crossed; in stop mode
        it stays where it was.
        """
        blocked, tangents = self.blocked_with_tangents(before, updates, pop)
        if len(blocked) == 0:
            return blocked
        self.total_blocked += len(blocked)
        start = before.loc[blocked, ["x", "y"]].to_numpy(dtype=float)
        if self.mode in ("slide", "reflect"):
            move = updates.loc[blocked, ["x", "y"]].to_numpy(dtype=float) - start
            along = (move * tangents).sum(axis=1, keepdims=True) * tangents
            velocity = updates.loc[blocked, ["vx", "vy"]].to_numpy(dtype=float)
            velocity_along = (velocity * tangents).sum(axis=1, keepdims=True) * tangents
            if self.mode == "reflect":
                # mirror the normal component: along - (move - along)
                updates.loc[blocked, ["x", "y"]] = start + 2 * along - move
                updates.loc[blocked, ["vx", "vy"]] = 2 * velocity_along - velocity
            else:
                updates.loc[blocked, ["x", "y"]] = start + along
                updates.loc[blocked, ["vx", "vy"]] = velocity_along
        else:
            updates.loc[blocked, ["x", "y", "z"]] = before.loc[
                blocked, ["x", "y", "z"]
            ].to_numpy()
        return blocked

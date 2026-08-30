"""Headless construction and execution of model specifications.

The model specifications in this repository include the pygame visualizer,
which is unhelpful for batch V&V runs. This module builds an
:class:`~vivarium.InteractiveContext` from a model specification YAML while
skipping visualization components.
"""

import importlib
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml
from vivarium import InteractiveContext

SKIPPED_MODULES = ("visualizer",)

TREE_ATTRIBUTES = [
    "x",
    "y",
    "z",
    "vx",
    "vy",
    "vz",
    "frozen",
    "path_id",
    "parent_id",
    "depth",
    "radius",
]


def build_headless_simulation(spec_path: str | Path) -> InteractiveContext:
    """Build an InteractiveContext from a model spec, skipping the visualizer."""
    with open(spec_path) as f:
        spec = yaml.safe_load(f)

    components = []
    for package, entries in spec.get("components", {}).items():
        for entry in entries:
            entry = entry.strip()
            if not entry.endswith("()"):
                raise NotImplementedError(
                    f"Only no-argument component constructors are supported, got: {entry}"
                )
            module_name, class_name = entry[:-2].rsplit(".", 1)
            if module_name.split(".")[0] in SKIPPED_MODULES:
                continue
            module = importlib.import_module(f"{package}.{module_name}")
            components.append(getattr(module, class_name)())

    configuration = spec.get("configuration", {})
    return InteractiveContext(components=components, configuration=configuration)


def get_ellipsoid_bounds(sim: InteractiveContext) -> tuple[float, float]:
    """Get the (a, b) semi-axes of the containment ellipsoid, if configured."""
    a, b, _ = get_ellipsoid_semi_axes(sim)
    return a, b


def get_ellipsoid_semi_axes(sim: InteractiveContext) -> tuple[float, float, float]:
    """Get the (a, b, c) semi-axes of the containment ellipsoid, if configured."""
    try:
        config = sim.configuration.ellipsoid_containment
        return float(config.a), float(config.b), float(config.c)
    except AttributeError:
        return 2.0, 2.0, 2.0


def get_network(sim: InteractiveContext) -> pd.DataFrame:
    """Get the particle table attributes needed for network analysis."""
    return sim.get_population(TREE_ATTRIBUTES)


def tree_edges(pop: pd.DataFrame) -> pd.DataFrame:
    """Get vessel segments as a frame of parent/child coordinate pairs.

    Segments connect each particle with a valid parent to that parent's
    position; particles whose parent has left the table are skipped.
    """
    children = pop[pop.parent_id >= 0]
    valid = children.parent_id.isin(pop.index)
    children = children[valid]
    parents = pop.loc[children.parent_id]
    return pd.DataFrame(
        {
            "x0": parents.x.values,
            "y0": parents.y.values,
            "z0": parents.z.values,
            "x1": children.x.values,
            "y1": children.y.values,
            "z1": children.z.values,
            "child": children.index.values,
            "parent": children.parent_id.values,
            "radius": children.radius.values,
        }
    )


def run_steps(sim: InteractiveContext, n_steps: int, callback=None, every: int = 1) -> None:
    """Step the simulation, optionally invoking a callback every ``every`` steps."""
    for step in range(n_steps):
        sim.step()
        if callback is not None and (step + 1) % every == 0:
            callback(step + 1)

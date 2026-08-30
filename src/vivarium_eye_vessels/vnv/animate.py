"""Render an animated GIF of the vessel formation process.

Usage::

    vnv_growth_gif src/vivarium_eye_vessels/model_specifications/model_spec.yaml

Runs the model specification headless (no pygame window) and draws one frame
every few time steps. By default the GIF overwrites ``docs/vnv/growth.gif``
in place, so committing it makes the before/after of a model change visible
as an image diff in the pull request.
"""

from pathlib import Path

import click
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from PIL import Image

from vivarium_eye_vessels.vnv import simulation

VESSEL_COLOR = "#2a78d6"  # categorical slot 1: the simulation's identity color
TIP_COLOR = "#16324f"
INK = "#333333"
MUTED = "#767676"


def render_frame(
    pop, bounds: tuple[float, float], day: float, size: int = 640
) -> Image.Image:
    """Draw the current vessel network as one animation frame."""
    a, b = bounds
    fig, ax = plt.subplots(figsize=(size / 100, size / 100), dpi=100)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    boundary = matplotlib.patches.Ellipse(
        (0, 0), 2 * a, 2 * b, fill=False, color="#dddddd", linewidth=1
    )
    ax.add_patch(boundary)

    edges = simulation.tree_edges(pop)
    if not edges.empty:
        segments = np.stack([edges[["x0", "y0"]].values, edges[["x1", "y1"]].values], axis=1)
        ax.add_collection(
            LineCollection(segments, colors=VESSEL_COLOR, linewidths=1.0, alpha=0.9)
        )

    tips = pop[~pop.frozen & (pop.path_id >= 0)]
    if not tips.empty:
        ax.scatter(tips.x, tips.y, s=6, color=TIP_COLOR, zorder=3)

    margin = 0.08
    ax.set_xlim(-a * (1 + margin), a * (1 + margin))
    ax.set_ylim(-b * (1 + margin), b * (1 + margin))
    ax.set_aspect("equal")
    ax.axis("off")

    n_frozen = int(pop.frozen.sum())
    ax.text(
        0.02,
        0.98,
        f"day {day:.1f}",
        transform=ax.transAxes,
        va="top",
        color=INK,
        fontsize=11,
    )
    ax.text(
        0.02,
        0.02,
        f"{len(pop)} particles · {n_frozen} vessel segments · {len(tips)} growth tips",
        transform=ax.transAxes,
        va="bottom",
        color=MUTED,
        fontsize=9,
    )

    fig.tight_layout(pad=0.3)
    fig.canvas.draw()
    frame = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert("RGB")
    plt.close(fig)
    return frame


@click.command()
@click.argument("model_spec", type=click.Path(exists=True))
@click.option(
    "-o", "--output", default="docs/vnv/growth.gif", show_default=True, type=click.Path()
)
@click.option("--steps", default=800, show_default=True, help="Total simulation steps.")
@click.option(
    "--steps-per-frame", default=8, show_default=True, help="Time steps per GIF frame."
)
@click.option("--fps", default=12, show_default=True, help="GIF frames per second.")
@click.option("--size", default=640, show_default=True, help="Frame size in pixels.")
def main(
    model_spec: str, output: str, steps: int, steps_per_frame: int, fps: int, size: int
) -> None:
    """Animate vessel formation for MODEL_SPEC as a GIF."""
    sim = simulation.build_headless_simulation(model_spec)
    bounds = simulation.get_ellipsoid_bounds(sim)
    step_size_days = float(sim.configuration.time.step_size)

    frames = [render_frame(simulation.get_network(sim), bounds, 0.0, size)]

    def capture(step: int) -> None:
        pop = simulation.get_network(sim)
        frames.append(render_frame(pop, bounds, step * step_size_days, size))

    simulation.run_steps(sim, steps, callback=capture, every=steps_per_frame)

    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Quantize for a small file; hold the final frame for a beat
    quantized = [f.quantize(colors=64, dither=Image.Dither.NONE) for f in frames]
    durations = [1000 // fps] * len(quantized)
    durations[-1] = 2000
    quantized[0].save(
        out_path,
        save_all=True,
        append_images=quantized[1:],
        duration=durations,
        loop=0,
        optimize=True,
    )
    click.echo(f"Wrote {len(frames)} frames to {out_path}")


if __name__ == "__main__":
    main()

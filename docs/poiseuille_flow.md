# Poiseuille flow, for people who model health rather than pipes

This is the physics background for the `FlowRemodeler` component (realism
roadmap idea 5). No fluid dynamics is assumed; the payoff is that one
19th-century formula turns our frozen particle table into a solvable
electrical circuit, and one biological observation turns the solution into
a rule for which vessels thicken and which ones die.

## One pipe

Push fluid through a thin straight tube slowly enough that it flows in
smooth layers (no turbulence — true for essentially all blood vessels
outside the largest arteries). The fluid at the wall doesn't move; the
fluid at the center moves fastest; the velocity profile in between is a
parabola. Working out the total volume passing through per second gives
the **Hagen–Poiseuille law** (Poiseuille was a French physician, 1840s,
studying exactly our problem — blood flow):

```
Q = π r⁴ Δp / (8 μ L)
```

- `Q` — volumetric flow (volume per second)
- `r` — the tube's inner radius
- `Δp` — the pressure drop from one end to the other
- `μ` — the fluid's viscosity ("thickness"; honey is high, water is low)
- `L` — the tube's length

The number to stare at is the **fourth power of the radius**. Double a
vessel's caliber and, for the same pressure drop, it carries *sixteen
times* the flow. Narrow a coronary artery by 20% and its capacity falls
by nearly 60%. This is why caliber is the quantity the body actively
regulates, why Murray's law is about radii, and why our `radius` column
is the most consequential one in the particle table.

## Ohm's law in disguise

Rearrange the formula as `Q = g · Δp` with

```
g = π r⁴ / (8 μ L)        (the "conductance" of the segment)
```

and it is exactly Ohm's law `I = V / R`: flow is current, pressure is
voltage, and every vessel segment is a resistor whose conductance you can
read off its radius and length. A vascular *network* is then just a
resistor network:

- **Kirchhoff's current law** at every junction: flow in = flow out
  (blood is incompressible and doesn't pile up at bifurcations).
- **Fixed voltages at the terminals**: arteries enter the retina at high
  pressure, veins leave at low pressure.

Writing Kirchhoff's law at every interior node gives one linear equation
per node — a sparse symmetric system (the weighted graph Laplacian, for
graph people). Solving it yields the pressure at every junction, and from
pressures you get the flow through every segment. At our scale (~10⁴
nodes) this is a millisecond-class sparse solve, which is why the roadmap
called it "cheap at this scale."

This is also why anastomosis (roadmap idea 4) had to come first: a tree
with no loops has no artery-to-vein path, so no current can flow. Circuits
perfuse; trees don't. In the simulation we also give every node a small
"leak" conductance to a virtual tissue node at intermediate pressure —
physically, the capillary drainage that happens everywhere along real
microvessels — so dead-end twigs carry a trickle rather than exactly zero,
and the linear system is always well-posed.

## Wall shear stress: what the vessel itself can feel

A vessel has no way to measure its own flow in ml/s. What its endothelial
lining *can* feel is the drag of blood rubbing past it — the **wall shear
stress**. For Poiseuille flow it is

```
τ = 4 μ Q / (π r³)
```

Two facts make this the master variable of vascular biology:

1. **Endothelial cells are shear sensors.** Sustained high shear makes a
   vessel remodel outward (widen); sustained low shear makes it regress
   and, below a threshold, get pruned away entirely. This is the
   Pries–Secomb structural adaptation model, and it is how development
   sculpts a functioning network out of an over-produced capillary mesh —
   the same over-produce-then-prune logic as synaptic pruning in the brain.

2. **Murray's law falls out for free.** If every vessel adapts until it
   feels the *same* target shear, then `Q ∝ r³` everywhere, and
   conservation of flow at a bifurcation (`Q₀ = Q₁ + Q₂`) becomes
   `r₀³ = r₁³ + r₂³` — exactly the junction exponent k = 3 that our V&V
   fits at the network's bifurcations. The splitter *imposes* Murray's law
   at birth; shear-driven remodeling is the mechanism that would *maintain*
   it in life.

## What the simulation does with all this

Every `remodel_interval` steps, `FlowRemodeler`:

1. builds the resistor network from the frozen particle table — one
   resistor per parent-child segment plus one per anastomosis bridge,
   conductance `r⁴ / L`, plus the per-node tissue leak;
2. fixes artery roots at `artery_pressure` and vein roots at
   `vein_pressure`, and solves the Kirchhoff system for all node pressures
   (sparse direct solve);
3. computes each segment's flow `Q = g Δp` and shear `τ ∝ Q / r³`;
4. **prunes**: terminal segments (degree-1 ends of the graph, so the
   network is never cut in the middle) whose shear is below
   `shear_threshold_fraction` × their own tree's median shear are
   recycled — the particle's `unfreeze_time` is stamped (the column was
   born for this) and it returns to the free pool as a fresh wandering
   particle;
5. **adapts** (optional): every segment's radius drifts by
   `adaptation_rate` toward the caliber that would put its shear at its
   own tree's median — high-shear trunks thicken, low-shear twigs thin.
   The target is per tree, not global, because arteries really do run at
   higher wall shear than veins (roughly 40–70 vs. 1–20 dyn/cm² in
   humans); one shared target would thicken the arteries until the
   clinical artery:vein caliber ratio inverted, which is exactly what
   happened on the first attempt.

Repeated passes eat dead-end branches back from their tips, which is
precisely the surgical tool for our biggest measured gap: the
diameter-stratified V&V showed the sim has ~20× too many capillary-scale
skeleton branches per image while matching real segment lengths within
every caliber stratum. Low-shear twig pruning removes exactly that clutter
while the loops idea 4 built — which carry real artery-to-vein flow —
are protected by their high shear.

## What we deliberately ignore (and why it's fine)

- **Blood is not a Newtonian fluid** — viscosity drops in small vessels as
  red cells single-file (the Fåhræus–Lindqvist effect). We use a constant
  μ, which cancels out of every *relative* comparison we make (shear
  vs. the network median).
- **Flow is pulsatile**, not steady. Microvascular remodeling responds to
  time-averaged shear, which the steady solution approximates.
- **Pressures are in arbitrary units** — only ratios matter for pruning
  decisions, so the boundary values are model knobs, not mmHg.

The classic references are Murray (1926) for the optimality argument and
Pries & Secomb (1998) for the adaptation model; both are in the roadmap's
reference list.

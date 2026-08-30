"""Validation and verification tools for the eye vessel simulation.

This subpackage provides:

- :mod:`~vivarium_eye_vessels.vnv.simulation`: run model specifications
  headless (without the pygame visualizer)
- :mod:`~vivarium_eye_vessels.vnv.metrics`: quantitative network metrics
  (fractal dimension, skeleton density, segment lengths, tortuosity,
  bifurcation angles)
- :mod:`~vivarium_eye_vessels.vnv.reference_data`: download and cache
  expert-labeled vessel masks from the public HRF dataset
- :mod:`~vivarium_eye_vessels.vnv.animate`: the ``vnv_growth_gif`` command
- :mod:`~vivarium_eye_vessels.vnv.compare`: the ``vnv_compare`` command

See ``docs/realism_roadmap.md`` for how these fit into the model development
workflow.
"""

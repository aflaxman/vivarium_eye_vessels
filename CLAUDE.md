# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Vivarium simulation model that creates synthetic data resembling the vascular system of the human eye. The simulation uses agent-based modeling with custom components for vessel growth, splitting, and collision avoidance. Components work together to simulate blood vessel development in 3D space using particle systems with physics-based interactions.

## Development Commands

### Testing
```bash
pytest                    # Run all tests
pytest tests/             # Run tests in tests directory  
pytest tests/test_sample.py  # Run specific test file
```

### Code Formatting and Linting
```bash
black .                   # Format code (line length: 94)
isort .                   # Sort imports (black profile)
black --check .           # Check formatting without changes
isort --check .           # Check import sorting without changes
```

### Running Simulations
```bash
simulate run src/vivarium_eye_vessels/model_specifications/model_spec.yaml
simulate run -v src/vivarium_eye_vessels/model_specifications/model_spec.yaml  # Verbose logging
```

### Artifact Generation
This sim does not use artifacts, actually

## Architecture

### Core Components (src/vivarium_eye_vessels/components/)

- **particles.py**: Core particle system with Particle3D base component, PathFreezer, PathSplitter, PathExtinction, and PathDLA for vessel growth dynamics
- **boundaries.py**: Force-based boundary conditions including EllipsoidContainment, CylinderExclusion (the caliber-aware foveal avascular zone), PointRepulsion, and FrozenRepulsion using Hookean or magnetic force calculations; growth guidance (PerfusionDemand hypoxia attraction, DevelopmentalWave growth front, ArcadeGuidance radial template for arcade-caliber tips)
- **visualizer.py**: 3D visualization using pygame for real-time particle rendering
- **observers.py**: Data collection and output management

### Key Vivarium Patterns (vivarium 4.x)

1. **Component Architecture**: All components inherit from `vivarium.Component` with standard lifecycle methods (`setup()`, `on_time_step()`, etc.)

2. **Configuration Structure**: Model specifications in YAML format define component parameters, with `CONFIGURATION_DEFAULTS` in each component class

3. **Population Management**: Particles are treated as simulants with tabular data (position, velocity, frozen state, parent relationships). `Particle3D` registers all particle columns via `builder.population.register_initializer(...)` and initializes them with `population_view.initialize(...)`

4. **Private Columns**: vivarium 4 only allows the component that created a column to write it. `Particle3D` owns all particle state and exposes `update_particles(updates)`; sibling components (PathFreezer, PathSplitter, PathExtinction, PathDLA) obtain it via `builder.components.get_components_by_type(Particle3D)[0]` and route their writes through it

5. **Explicit Attribute Reads**: `population_view.get(index, attributes)` requires an explicit attribute list; components declare theirs in a `required_attributes` property

6. **Event-Driven Updates**: Components respond to time step events to update particle states and apply forces

7. **Builder Pattern**: Use `Builder` object in `setup()` methods to register value pipelines (`register_value_producer`/`register_value_modifier` with `required_resources`), event listeners, and simulant initializers

### Data Structure

- **Particle Columns**: x, y, z (position), vx, vy, vz (velocity), frozen, freeze_time, unfreeze_time, depth, parent_id, path_id
- **Force Types**: Hookean (spring-based) and magnetic (inverse square) force calculations
- **Tree Structure**: Vessels maintain parent-child relationships for branching patterns

### Configuration Management

Model specifications use nested YAML with component-specific parameter sections. Critical parameters include:
- Particle dynamics (velocity limits, force thresholds)  
- Boundary constraints (ellipsoid/cylinder geometry, spring constants)
- Vessel growth (split intervals, angles, probabilities)
- Visualization settings (colors, projection scale, frame rate)

### Code Style

- Python 3.10/3.11 support
- Type hints required for all functions and methods
- Black formatting (94 character line length)  
- Sparse comments only for complex operations
- Descriptive variable and function names to minimize comment needs
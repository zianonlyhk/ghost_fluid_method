# Configuration Examples

This directory contains YAML configuration files for different types of Ghost Fluid Method simulations. Each example demonstrates specific flow phenomena and numerical techniques.

## Configuration Format

All examples use YAML format (.yaml) which provides:
- **Syntax highlighting** in modern editors (VS Code, Vim, Sublime, etc.)
- **Better readability** with structured layout
- **Error detection** and validation support
- **IntelliSense/autocomplete** in many editors

## Available Examples

### 1. Shock-Circle Interaction (`shock_circle_interaction.yaml`)
**Phenomenon**: High-pressure shock wave hitting a stationary circular obstacle
- **Physics**: Shock diffraction, reflection, vortex formation
- **Grid**: 400×200 cells for shock resolution
- **Key Features**: Strong density/pressure ratios, fixed circular obstacle
- **Runtime**: ~5-10 minutes on modern hardware

### 2. Flow Past Cylinder (`flow_past_cylinder.yaml`)
**Phenomenon**: Uniform subsonic flow encountering a stationary cylinder
- **Physics**: Flow separation, wake formation, potential vortex shedding
- **Grid**: 300×150 cells for wake capture
- **Key Features**: Subsonic flow, long domain for wake development
- **Runtime**: ~8-15 minutes on modern hardware

### 3. Moving Body in Fluid (`moving_body_fluid.yaml`)
**Phenomenon**: Accelerating rigid body through initially quiescent fluid
- **Physics**: Fluid-structure interaction, added mass effects, body-generated flow
- **Grid**: 200×200 cells (square domain)
- **Key Features**: Enabled rigid body acceleration, initially stationary fluid
- **Runtime**: ~6-12 minutes on modern hardware

### 4. Double Circle Complex (`double_circle_complex.yaml`)
**Phenomenon**: Flow through complex double-circle geometry
- **Physics**: Flow channeling, multiple obstacle interactions, geometric complexity
- **Grid**: 300×225 cells
- **Key Features**: Complex geometry definition, moderate flow conditions
- **Runtime**: ~10-18 minutes on modern hardware

## How to Use Examples

### Method 1: Direct Execution
```bash
# Run example directly (no copying needed)
./bin/main examples/shock_circle_interaction.yaml
./bin/main examples/flow_past_cylinder.yaml
./bin/main examples/moving_body_fluid.yaml
```

### Method 2: Copy and Customize
```bash
# Copy example configuration and modify
cp examples/shock_circle_interaction.yaml my_config.yaml
code my_config.yaml  # Edit with syntax highlighting
./bin/main my_config.yaml
```

### Method 3: Edit in Place
```bash
# Edit example directly with syntax highlighting
code examples/shock_circle_interaction.yaml
vim examples/flow_past_cylinder.yaml
# Then run directly
./bin/main examples/shock_circle_interaction.yaml
```

### Method 4: Parameter Studies
```bash
# Create a study directory
mkdir parameter_study
cp examples/flow_past_cylinder.yaml parameter_study/base_config.yaml

# Modify parameters for different runs with YAML syntax
# Example: vary Reynolds number by changing density/velocity
./bin/main parameter_study/base_config.yaml
```

## Customization Guidelines

### Grid Resolution
- **Low resolution** (50×50 to 100×100): Quick testing, parameter exploration
- **Medium resolution** (150×150 to 300×300): Production runs, reasonable accuracy
- **High resolution** (400×400+): High-accuracy studies, publication quality

### CFL Number Selection
- **c = 0.5-0.7**: Very stable, slower convergence
- **c = 0.8-0.9**: Standard choice, good balance
- **c = 0.95-1.0**: Aggressive, may cause instability

### Time Duration Guidelines
- **Short runs** (tStop = 0.5-1.0): Initial testing, debugging
- **Medium runs** (tStop = 2.0-5.0): Most physical phenomena
- **Long runs** (tStop = 10.0+): Steady state, statistical averaging

### Output Frequency
- **loggingFactor = 1-5**: Detailed time evolution, large files
- **loggingFactor = 10-20**: Standard temporal resolution
- **loggingFactor = 50+**: Minimal output, final state focus

## Performance Estimates

Runtime estimates on modern desktop (Intel i7/AMD Ryzen 7):

| Grid Size | Time Steps | Estimated Runtime |
|-----------|------------|-------------------|
| 100×100   | 1000       | 1-2 minutes      |
| 200×200   | 2000       | 5-10 minutes     |
| 400×400   | 4000       | 30-60 minutes    |
| 600×600   | 6000       | 2-4 hours        |

**Memory Usage**: Approximately 50-100 MB per 100×100 grid

## Validation and Reference Data

Each example can be compared against:
- **Literature results**: See main README references
- **Analytical solutions**: Available for simple geometries
- **Commercial CFD**: ANSYS Fluent, OpenFOAM comparisons
- **Experimental data**: Wind tunnel, shock tube experiments

## Troubleshooting Examples

### Common Issues
1. **Simulation crashes**: Reduce CFL number, check boundary conditions
2. **Unphysical results**: Verify initial conditions, check rigid body placement
3. **Slow convergence**: Increase grid resolution, adjust reinitialization
4. **Large file sizes**: Increase loggingFactor, reduce simulation time

### Quick Debugging
```bash
# Test with minimal resolution first (YAML format)
sed -i 's/nCells_x: [0-9]*/nCells_x: 50/' config.yaml
sed -i 's/nCells_y: [0-9]*/nCells_y: 50/' config.yaml
sed -i 's/tStop: [0-9.]*/tStop: 0.1/' config.yaml
```

## Creating New Examples

When creating custom configurations:

1. **Start from similar example**: Copy closest matching case
2. **Modify incrementally**: Change one parameter at a time
3. **Test at low resolution**: Verify setup before high-resolution runs
4. **Document changes**: Add comments explaining modifications
5. **Validate results**: Compare with known solutions when possible
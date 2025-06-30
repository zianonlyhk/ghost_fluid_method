# Ghost Fluid Method for 2D Euler Equations

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## About

This project implements the Ghost Fluid Method (GFM) for solving 2D Euler equations with fluid-structure interaction. The code is written in C++ and uses a finite volume approach with MUSCL-Hancock HLLC scheme for high-resolution shock capturing.

**Key Features:**
- Ghost Fluid Method for immersed solid boundaries
- MUSCL-Hancock HLLC Riemann solver
- Level set method for interface tracking
- Rigid body motion with fluid coupling
- Multiple geometric shapes (circle, square, double circle)
- Configurable initial conditions and boundary states

## Getting Started

### Prerequisites
- **GNU Make**: Build system
- **C++ Compiler**: g++ or clang++ with C++17 support
- **Gnuplot**: For visualization and animation generation

### Building the Project

1. **Clone the repository:**
```bash
git clone https://github.com/[username]/ghost_fluid_method.git
cd ghost_fluid_method
```

2. **Compile the code:**
```bash
make clean && make
```

This will build the executable `bin/main` from source files in `src/`.

**Build Options:**
- `make` or `make release`: Optimized release build (-O3)
- `make clean`: Remove object files and executable
- `make clean_all`: Remove all generated files including output data

### Running the Solver

**Basic Usage:**
```bash
./bin/main [config_file.yaml]
```

**Examples:**
```bash
./bin/main                                    # Use default config.yaml
./bin/main examples/shock_circle_interaction.yaml  # Use example configuration
./bin/main my_simulation.yaml                # Use custom configuration
./bin/main --help                            # Show help message
```

## Configuration Guide

The Ghost Fluid Method solver uses YAML format for configuration files, providing excellent syntax highlighting and readability in modern text editors.

### Configuration File Format

**YAML Format:**
```yaml
# Comments start with #
parameter_name: value
array_parameter: [value1, value2, value3]
```

**Benefits of YAML:**
- ✅ **Syntax highlighting** in all modern editors (VS Code, Vim, Sublime, etc.)
- ✅ **Better readability** with structured layout
- ✅ **Error detection** and validation support
- ✅ **IntelliSense/autocomplete** in many editors
- ✅ **Cleaner array syntax**: `[1.0, 2.0, 3.0]`
- ✅ **Type safety** with automatic validation

### Domain Configuration

**Grid Boundaries:**
```yaml
x0: 0.0         # Left boundary x-coordinate
x1: 1.0         # Right boundary x-coordinate  
y0: 0.0         # Bottom boundary y-coordinate
y1: 1.0         # Top boundary y-coordinate
```
- Defines the physical domain [x0,x1] × [y0,y1]
- Must satisfy: x0 < x1 and y0 < y1
- Typical values: Unit square [0,1] × [0,1]

**Grid Resolution**
```
nCells_x 241    # Number of cells in x-direction
nCells_y 241    # Number of cells in y-direction
```
- Higher resolution = more accurate but slower simulation
- Odd numbers often work better for centered objects
- Typical range: 50-500 cells per direction
- Memory usage scales as nCells_x × nCells_y

### Time Integration Settings

**CFL Number**
```
c 0.9           # Courant-Friedrichs-Lewy number
```
- Controls time step size for numerical stability
- Must be: 0 < c ≤ 1.0
- Lower values = more stable but slower
- Recommended: 0.8-0.9 for most cases

**Simulation Duration**
```
tStop 2.0       # Simulation end time
```
- Physical time when simulation stops
- Units depend on your problem scaling
- Typical values: 0.1-10.0

### Rigid Body Configuration

**Position and Velocity**
```
initRigidBodyLoc 0.2 0.5     # Initial (x,y) center coordinates
initRigidBodyVel 0.0 -1.0    # Initial (vx,vy) velocity
```
- Position must be within domain boundaries
- Velocity in physical units

**Rigid Body Type**
```
rigidBodyType 1              # 1=circle, 2=square, 3=double circle
rigidBodyRadiusOrLength 0.125    # Size parameter
rigidBodyAdditionalFactor 0.0    # Extra parameter (spacing for double circles)
```

**Rigid Body Types:**
- **Type 1 (Circle)**: `rigidBodyRadiusOrLength` = radius
- **Type 2 (Square)**: `rigidBodyRadiusOrLength` = side length  
- **Type 3 (Double Circle)**: Two circles with separation = `rigidBodyAdditionalFactor`

**Motion Control**
```
acceleratingRigidBody 1      # 1=enable acceleration, 0=disable
```
- Enables/disables rigid body acceleration due to fluid forces
- Useful for studying fluid-structure interaction

### Initial Fluid Conditions

**Left and Right States**
```
leftState 1.0 0.0 0.0 1.0    # (density, u-velocity, v-velocity, pressure)
rightState 1.0 0.0 0.0 1.0   # (density, u-velocity, v-velocity, pressure)
leftRightStateBoundary 0.2   # x-coordinate separating left/right states
```

**State Vector Format: [ρ, u, v, p]**
- **ρ (density)**: Mass per unit volume (typically 0.1-10.0)
- **u (x-velocity)**: Horizontal velocity component
- **v (y-velocity)**: Vertical velocity component  
- **p (pressure)**: Fluid pressure (typically 0.1-10.0)

**Common Initial Conditions:**
- **Shock tube**: High pressure left, low pressure right
- **Uniform flow**: Same state on both sides with non-zero velocity
- **At rest**: Zero velocities, equal pressures

### Numerical Parameters

**Level Set Reinitialization**
```
reinitFactor 1               # Reinitialize level set every N steps
```
- Maintains level set as signed distance function
- Lower values = more accurate but slower
- Typical range: 1-10

### Output Configuration

**Logging Frequency**
```
loggingFactor 15             # Write output every N time steps
```
- Controls how often results are saved
- Lower values = more output files, larger storage
- Balance between temporal resolution and file size

**Output Naming**
```
runName gfm_sim              # Base name for output files
```
- Prefix for all output file names
- Use descriptive names for different test cases

## Running Simulations

### Quick Start

1. **Configure your simulation**:
   ```bash
   # Edit YAML configuration with syntax highlighting
   code config.yaml
   # Or use any YAML-capable editor
   vim config.yaml
   ```

2. **Build the solver**:
   ```bash
   make clean && make
   ```

3. **Run the simulation**:
   ```bash
   ./bin/main                                   # Use default config.yaml
   # OR specify a custom configuration file
   ./bin/main examples/shock_circle_interaction.yaml
   ./bin/main my_custom_config.yaml
   ```

### Simulation Workflow

The solver executes the following steps:

1. **Configuration Loading**: Reads and validates `config.yaml`
2. **Domain Initialization**: Sets up computational grid and ghost cells
3. **Initial Conditions**: Applies fluid states and rigid body placement
4. **Main Loop**: Advances solution in time using:
   - Ghost fluid method for solid boundaries
   - MUSCL-Hancock HLLC scheme for fluid dynamics
   - Level set advection and reinitialization
   - Rigid body motion (if enabled)
5. **Output Generation**: Writes results to `./output/`

### Monitoring Progress

During execution, the solver displays:
```
1: 0 / 2.0
2: 0.0123 / 2.0
3: 0.0247 / 2.0
...
```
Format: `[output_step]: [current_time] / [end_time]`

### Output Files

Results are saved in `./output/` with the following naming convention:
- `{runName}_rhoResults.dat` - Density field
- `{runName}_pressureResults.dat` - Pressure field  
- `{runName}_velMagResults.dat` - Velocity magnitude
- `{runName}_levelSetResults.dat` - Level set function
- `{runName}_itnEnergyResults.dat` - Internal energy
- `{runName}_msResults.dat` - Mock Schlieren

### Common Test Cases

**1. Shock-Circle Interaction:**
```yaml
# Moving shock hits stationary circle
leftState: [1.4, 0.0, 0.0, 1.4]    # High pressure region
rightState: [1.0, 0.0, 0.0, 1.0]   # Low pressure region
leftRightStateBoundary: 0.2
rigidBodyType: 1
rigidBodyRadiusOrLength: 0.1
acceleratingRigidBody: 0            # Fixed circle
```

**2. Flow Past Cylinder:**
```yaml
# Uniform flow past cylinder  
leftState: [1.0, 0.5, 0.0, 1.0]    # Moving fluid
rightState: [1.0, 0.5, 0.0, 1.0]   # Same state
leftRightStateBoundary: 0.1         # Small left region
rigidBodyType: 1
acceleratingRigidBody: 0            # Fixed cylinder
```

**3. Moving Body in Fluid:**
```yaml
# Accelerating body through stationary fluid
leftState: [1.0, 0.0, 0.0, 1.0]    # Stationary fluid
rightState: [1.0, 0.0, 0.0, 1.0]
initRigidBodyVel: [0.0, -0.5]      # Initial downward motion
acceleratingRigidBody: 1            # Allow acceleration
```

## Code Structure

### Main Components

**Core Solver**
- `src/main.cc`: Main driver program with simulation loop
- `src/gfm_2d_euler_solver.cc/hh`: Core Ghost Fluid Method implementation
- `src/flux_func.cc/hh`: Riemann solver and flux calculations
- `src/ghost_fluid_utilities.cc/hh`: GFM-specific boundary operations

**Configuration and Setup**
- `src/config_manager.cc/hh`: Configuration file parsing and validation
- `src/initial_conditions.cc/hh`: Initial condition setup
- `src/level_set_functions.cc/hh`: Level set geometry functions
- `src/constants.hh`: Physical and numerical constants

**Utilities**
- `src/vec_transform.cc/hh`: Vector/matrix operations
- `src/inline/`: Inline utility functions
  - `primitive_tran.hh`: Primitive/conservative variable transformations
  - `cell_operation.hh`: Grid cell operations
  - `debug_tools.hh`: Debugging and visualization helpers

### Build System
- `Makefile`: Automated build configuration
- Supports debug and release builds
- Automatic dependency tracking

## Output Format and Visualization

### Data File Structure

Simulation results use a structured ASCII format optimized for Gnuplot:

```
time x y value
```

**Format Rules:**
- Each line contains one data point: `time x_coordinate y_coordinate field_value`
- Single blank line separates different y-coordinates (new row)
- Double blank line separates time steps
- Data is organized in row-major order (x varies fastest)

**Example Data Structure:**
```
# Time step 1 (t=0.0)
0.0 0.0 0.0 1.0   # Point (0,0)
0.0 0.1 0.0 1.0   # Point (0.1,0)
0.0 0.2 0.0 1.0   # Point (0.2,0)
                  # Blank line = next row
0.0 0.0 0.1 1.0   # Point (0,0.1)
0.0 0.1 0.1 1.0   # Point (0.1,0.1)
0.0 0.2 0.1 1.0   # Point (0.2,0.1)

                  # Double blank = next time step
# Time step 2 (t=0.05)
0.05 0.0 0.0 0.98
0.05 0.1 0.0 0.99
...
```

### Loading Data in Analysis Tools

**Python (with NumPy/Matplotlib):**
```python
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt('output/gfm_sim_pressureResults.dat')
t, x, y, p = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Plot final time step
final_time = np.max(t)
final_data = data[t == final_time]
plt.scatter(final_data[:, 1], final_data[:, 2], c=final_data[:, 3])
plt.colorbar(label='Pressure')
plt.show()
```

**MATLAB:**
```matlab
% Load and reshape data
data = load('output/gfm_sim_pressureResults.dat');
t = data(:,1); x = data(:,2); y = data(:,3); p = data(:,4);

% Extract final time step
final_idx = (t == max(t));
scatter(x(final_idx), y(final_idx), 50, p(final_idx), 'filled');
colorbar; title('Pressure Field');
```

### Visualization with Gnuplot

The included `plot_gif.gp` script creates animated visualizations:

```bash
# Generate animations after simulation
gnuplot plot_gif.gp
```

**Generated Outputs:**
1. Individual frame plots in `./plots/`
2. Animated GIFs in `./output/gif/`:
   - `{runName}_pressure.gif` - Pressure field evolution
   - `{runName}_rho.gif` - Density field evolution  
   - `{runName}_velMag.gif` - Velocity magnitude
   - `{runName}_mockschlieren.gif` - Mock Schlieren visualization

**Customizing Visualizations:**

Edit `plot_gif.gp` to modify:
- Color schemes and ranges
- Plot resolution and frame rate
- Field variables to visualize
- Animation speed and quality

### Advanced Analysis

**Extract Specific Data:**
```bash
# Get pressure along centerline at final time
awk '$1==max_time && $3==0.5 {print $2, $4}' output/gfm_sim_pressureResults.dat > centerline_pressure.dat
```

**Compare Multiple Runs:**
```bash
# Organize results by test case
mkdir results/shock_test_1 results/shock_test_2
mv output/*Results.dat results/shock_test_1/
# Run second case...
mv output/*Results.dat results/shock_test_2/
```

### Performance and Storage

**File Sizes:**
- Grid: 100×100, Time steps: 100 → ~40MB per field
- Grid: 200×200, Time steps: 200 → ~320MB per field

**Optimization Tips:**
- Increase `loggingFactor` to reduce output frequency
- Use shorter `tStop` for parameter studies  
- Clean old results: `make clean_all`

## Troubleshooting

### Common Configuration Errors

**"Missing required parameters"**
- Check that all required parameters are present in `config.yaml`
- Verify parameter names match exactly (case-sensitive)
- Ensure no typos in parameter names
- Use YAML syntax: `parameter: value` not `parameter value`

**"leftRightStateBoundary must be between x0 and x1"**
- Verify: `x0 < leftRightStateBoundary < x1`
- Check domain boundaries make sense

**"Value must be positive"**
- Grid sizes (`nCells_x`, `nCells_y`) must be > 0
- Time parameters (`tStop`, `loggingFactor`) must be > 0
- CFL number must be between 0 and 1

### Simulation Issues

**Simulation crashes/becomes unstable:**
- Reduce CFL number (try `c 0.5`)
- Check initial conditions for extreme values
- Verify rigid body is within domain
- Increase grid resolution

**Results look wrong:**
- Check initial conditions and boundary placement
- Verify rigid body type and size parameters
- Review visualization color scales
- Compare with reference solutions

**Slow performance:**
- Reduce grid resolution for testing
- Increase `loggingFactor` to output less frequently
- Use release build: `make clean && make`

### Getting Help

- Check configuration examples in `tmp/config_files/`
- Review test results in `tmp/results_configs/`
- Compare with reference papers and validation cases

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References

**Core Method:**
- Sambasivan, S. K., & UdayKumar, H. S. (2009). Ghost Fluid Method for strong shock interactions Part 2: Immersed Solid boundaries. AIAA Journal, 47(12), 2923–2937. https://doi.org/10.2514/1.43153

**Numerical Methods:**
- Toro, E. F. (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics. Springer.
- Osher, S., & Fedkiw, R. (2003). Level Set Methods and Dynamic Implicit Surfaces. Springer.

**Validation Cases:**
- Woodward, P., & Colella, P. (1984). The numerical simulation of two-dimensional fluid flow with strong shocks. Journal of Computational Physics, 54(1), 115-173.
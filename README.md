# Ghost Fluid Method for 2D Euler Equations

## About <a name = "about"></a>

This project implements the Ghost Fluid Method (GFM) for solving 2D Euler equations with fluid-structure interaction. The code is written in C++ and uses a finite volume approach.

## Getting Started <a name = "getting_started"></a>

### Prerequisites
- GNU Make
- C++ compiler (g++ or clang++)
- Gnuplot (for visualization)

### Building the Project

Compile the code using:
```
make
```

This will build the executable from source files in `src/`.

## Configuration <a name = "usage"></a>

Edit `config.txt` to configure simulation parameters. The key parameters include:
- Grid resolution
- Time step settings
- Physical constants
- Boundary conditions
- Output directory

## Running Simulations

After building, run the simulation with:
```
./bin/main
```

The solver will:
1. Read configuration from `config.txt`
2. Run the simulation
3. Output results to `./output`

## Code Structure

Key components:
- `src/main.cc`: Main driver program
- `src/gfm_2d_euler_solver.cc/hh`: Core solver implementation
- `src/flux_func.cc/hh`: Riemann solver and flux calculations
- `src/ghost_fluid_utilities.cc/hh`: GFM-specific operations
- `src/vec_transform.cc/hh`: Vector/matrix operations

## Visualization

After running the simulation with `./bin/main`, generate animations with:
```
gnuplot plot_gif.gp
```

This will:
1. Process simulation output from `./output`
2. Generate visualization frames in `./plots`
3. Create an animated GIF of the results

## Output Format

Simulation results are stored with the following format:
```
time x y value
```

Where:
- Blank lines separate y-coordinates
- Double blank lines separate time steps

Example for pressure at (x=0.2, y=0.5) at t=0.0:
```
0.0 0.2 0.5 1.0
```

A series of data arrays increments in the x-direction, with a newline indicates a jump alongs the y-axis. A double newline refers to a timestep forward, as shown in the following example:
```
...
0.0 0.3 0.4 0.9
0.0 0.4 0.4 0.9
0.0 0.5 0.4 0.9

0.0 0.0 0.5 0.9
0.0 0.1 0.5 0.9
0.0 0.2 0.5 0.9
0.0 0.3 0.5 0.9
0.0 0.4 0.5 0.9
0.0 0.5 0.5 0.9


0.0 0.0 0.0 1.0
0.0 0.1 0.0 1.0
...
```

## License <a name = "license"></a>

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References
- Sambasivan, S. K., & UdayKumar, H. S. (2009). Ghost Fluid Method for strong shock interactions Part 2: Immersed Solid boundaries. AIAA Journal, 47(12), 2923–2937. https://doi.org/10.2514/1.43153

# Fluid-structure interaction using the Ghost Fluid Method

## About <a name = "about"></a>

This project is part of the written assignment in the MPhil in Scientific Computing at the University of Cambridge, produced by student 2023P3. This file contains instructions on how one can reproduce the results presented in the final report.

## Unzipping <a name = "getting_started"></a>

Download from the repository to have the zip file "`MPhil_writtenAssignment_GFM-main.zip`".

Unzip by first changing to the directory containing the zip file and run:

```
unzip ./MPhil_writtenAssignment_GFM-main.zip
```

## Editting config file <a name = "usage"></a>

Going to the unzipped directory and edit the file "`./config.txt`":
```
cd ./MPhil_writtenAssignment_GFM-main/
vi ./config.txt
```

Edit the last line and change the "`/foo/bar`" part into the repository directory:
`````
repoDir /foo/bar/MPhil_writtenAssignment_GFM-main/
`````

For example, with the repository under the user home directory, this line is turned into:
`````
repoDir /home/2023P3/MPhil_writtenAssignment_GFM-main/
`````

## Compilation and execution <a name = "usage"></a>

Compile the binary executable by going to the unzipped directory and using the Make software:
```
make
```

The compiled binary is located under "`./bin/`". Start the simulation by running:
```
./bin/run_simulation
```

## Reproducing results

All initial conditions setting files are inside the "`./config_files/`" directory. To reproduce the simulation, simply copy the content in the file except the last line to the "`config.txt`" to overwrite configuration settings. 

When switching to the reflective boundary conditions, the source code file "`./src/exec_interface`" needs to be edited and the binary is required to be recompiled. The reflective boundary conditions can be achieved by commenting out the transmissive update functions and decommenting the refective functions inside the main while loop.

## Simulation data and visualisation

After running the simulation, the results are stored under the "`./data/`" directory. Different physical quantities are recorded in the format "`*Results_dat`".

Each line in a results file has the structure:
```
time x-coor y-coor quantity
```

For a cell at location `(x,y)=(0.2,0.5)`, having a pressure `p=1.0` at time `t=0.0`, the corresponding line in the "`*_pressureResults.dat`" is expressed as:
```
0.0 0.2 0.5 1.0
```

A blank line is used to indicate a jump to the next row along the y-axis. Two blank lines are used to locate a jump to the next time frame. These separation techniques can be picked up by the two features "`index`" and "`using 2:3:4 with image`" in gnuplot for visualising the simulation.

A collection of gnuplot scripts is included in the repository. These scripts are not prepared to run on all platforms. Extra work is requried to adapt these codes to one's use.
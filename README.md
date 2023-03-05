# Fluid-structure interaction using the Ghost Fluid Method

## About <a name = "about"></a>

This project is part of the written assignment in the MPhil in Scientific Computing at the University of Cambridge, by student 2023P3.

## Unzipping <a name = "getting_started"></a>

Download the repository to have the zip file "`MPhil_writtenAssignment_GFM-main.zip`".

Unzip by first changing to the directory containing the zip file and run:

```
unzip ./MPhil_writtenAssignment_GFM-main.zip
```

## Editting config file <a name = "usage"></a>

Going to the unzipped directory and edit the file "`./config.txt`":
```
cd MPhil_writtenAssignment_GFM-main
vi config.txt
```

Edit the last line and change the "`/foo/bar`" part into the repository directory:
`````
repoDir /foo/bar/MPhil_writtenAssignment_GFM-main/
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
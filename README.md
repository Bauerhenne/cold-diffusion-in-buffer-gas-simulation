# Cold-diffusion-in-buffer-gas-simulation
Simulation of the cooling of naphthalene molecules through a buffer gas in a spherical cell.

## Compilation

Compiling the library using make and a c++ compiler just with

```bash
make
```

The default c++ compiler is g++. If you want to use a different compiler, change the Makefile accordingly.

## Usage

In order to perform, one needs the file containing the field information of the buffergas on a discretized grid. The name of this file is given as first argument during calling the simulation program:

```bash
./simulation voxelized_D16mm_H8mm_d1.5mm_0deg_12Pa_50mW.txt
```

Above, the file is called "voxelized_D16mm_H8mm_d1.5mm_0deg_12Pa_50mW.txt". 

But the name can be arbitrarily chosen, but it should start with "voxelized_". 

The first line of the file starting with a "#" provides the diameter, the height, the exitDiameter, the delta (of the grid used) in the unit m, and the number of points in x,y,z direction. The following lines provide the velocity vx, vy, vz, (unit m/s) the density rho (unit particle/m³) and the temperature (unit K) of the buffergas on all grid points.

The simulation program calculates the trajectories of the naphthalene molecules in the buffergas within the cell and saves the information about each individual trajectory in the file naphthalene_step_NR.txt where NR is an ascending integer number. If the molecule hits the wall within a distance of 0.003m from the entry, the trajectory is not saved and the warning DISTANCE is to small! is printed by the program to stdout, where DISTANCE is the reached distance. The first two lines of the trajectory file start with "#" and provide information about the columns. Then each line provides the total number of collisions with the buffergas, the time, then the cartesian x,y,z coordinates, the x,y,z velocity of the naphthalene molecule. For the first 200 collisions, every timestep is printed, then only every 1000th timestep is printed.

In addition, the simulation program writes the file 

```bash
run_infos_${NAME}.txt
```

where NAME is the name of the voxelized file without the starting string "voxelized_".

This file contains in the first line starting with "#" a header describing the content of the columns. Then for each saved trajectory a line is added. The first number is the runNr, an ascending integer number. Then, if the trajectory ends within the buffergas sphere, a "0" is written, if it leaves the exit, then a "-1" is written. Then the number of timesteps and the time of flight of the trajectory is written, then the total number of collisions and then the x,y,z coordinate of the molecule at the end of the trajectory.

The program load_and_convert_u_data generates the voxelized file from a data file created from Comsol. This program takes the name of the Comsol file as argument. The filename must have the following format: 

```bash
D${NR1}mm_H${NR2}mm_d${NR3}mm_${something}.txt
```

The first number NR1 is the diameter, the second NR2 the height, the third NR3 the exit diameter. A valid file name is D16mm_H8mm_d1.5mm_0deg_12Pa_50mW.txt

The program writes the voxelized data in one output file with the name

````bash
voxelized_${filename}
````

where filename is the name of the Comsol file given as an input argument. 

## Data available

The zip file comsol.zip contains several data files containing buffer gas field information generated with Comsol. The zip file voxelized.zip contains the voxelized files generated by the program load_and_convert_u_data from the Comsol files in the comsol.zip file. These voxelized files can be directly used by the simulation program.
# Parallel Heat Transfer 2D
A classic MPI communication scheme for domain decomposition in 2D structured grids

## The problem

The 2D heat equation is defined as:
$$\frac{\delta u}{ \delta t} = \alpha (\frac{\delta u}{\delta x} + \frac{\delta u}{ \delta y})$$

when $\frac{\delta u}{\delta t}=0$ a steady-state heat transfer is attained at thermal equilibrium.
This simulation aims at attaining that steady state using a control volume discretization on a cartesian grid.

## Domain decomposition

The problem consists of a 2D domain, in which the left and right temperatures differ. The top and bottom of the plate are considered adiabatic, so no heat is transferred through them.

For parallelization, the domain is partitioned into subdomains that are assigned to the MPI processes.

![domain_decomposition](https://github.com/gkigiermo/parallel-heat-transfer/blob/main/imgs/domain_decomposition.png)

MPI Communications are used to maintain data coherence.

## Compilation
The compilation requires to have installed MPI in one of its flavors. 

To compile just use the makefile
```
make
```
that generates the executable ```par-heat-tranfer-2d.x```

## Execution
The execution must be done using MPI to enable the parallelization

```
mpirun -np 4 ./par-heat-tranfer-2d.x
```

The number of processes can change, but you will need to change the variables  ```npx, npx``` (```HeatTransfer2D.cpp```) to match the number of processes requested.
The correct execution displays the algorithm convergence and at the end generates the file ```temperatures.dat```

## Visualization
For visualizing the results you need to install Tecplot and load the generated file.
With the current configuration, the results should look like this:

<img src="https://github.com/gkigiermo/parallel-heat-transfer/blob/main/imgs/result.jpg" width="400">

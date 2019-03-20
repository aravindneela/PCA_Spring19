CMT-bone-BE is a skeleton app of CMT-nek, a Computational Fluid Dynamics (CFD) application. 
It is by default an MPI application that runs on a minimum of 2 process.
The computational fluid domain is partitioned in a 3D mesh fashion forming small cubical elements.
Many such elements are run on many processors using MPI where lots of computations are performed on each element. 
After the computation, MPI communication also takes place to exchange face information between adjacent elements.

Compiling the code:
$ make

Run the code:
$ mpirun -np <# of processors> ./cmtbonebe <timesteps> <Element size> <Element_x> <Element_y> <Element_z> <cart_x> <cart_y> <cart_z>

Here is a breif description of the command line arguments:
Timesteps: Number of times you want the simulation to run for. Run the simulation for 50 to 100 timesteps.
Element size: The size of each element in the 3D mesh. The range is 5 to 25. 
Element_x: Number of elements per processor along x-axis
Element_y: Number of elements per processor along y-axis
Element_z: Number of elements per processor along z-axis
cart_x: Number of processors along x-axis
cart_y: Number of processors along y-axis
cart_z: Number of processors along z-axis

The total number of elements per processor is Element_x * Element_y * Element_z
The total number of processors is cart_x * cart_y * cart_z.
NOTE: Make sure that <# of processors> (given through the -np flag in mpirun) is equal to cart_x * cart_y* cart_z.
      If this is not maintained, you will get a runtime error.

# Project5
This project is solved by creating algorithms for the explicit and implicit Euler scheme, as well as for the implicit Crank-Nicolson scheme. 
These algorithms are implemented in a C++ program, which is parallelized using OpenMP. This project also includes a program written in Python, which reads data from
output files, before plotting the results. The programs are executed as follows:

## One-dimensional diffusion equation, project5 directory
### main.cpp

The program takes two input arguments from the command line:
- 1: the number of timesteps
- 2: decides which scheme to use, varies between 1 and 4. 1 solves the explicit Euler scheme, 2 solves the implicit Euler scheme, 3 solves the implicit Crank-Nicolson scheme and 4 solves the explicit Euler scheme for unstable solutions.

## Two-dimensional diffusion equation, project5_mpi directory
### main.cpp
To compile this program on a Mac, write "clang-omp++ -fopenmp -03 main.cpp solvers.cpp -o project5_mp.x" in a terminal. 
To test for stability, please change the "dx"-value manually.

The program takes three input arguments from the command line:
- 1: the dimensions of the mesh
- 2: the number of timesteps
- 3: decides which solver to use, write "1" for Jacobi and "2" for Gauss-Seidel

## Python file, files_project5 directory
### project5.py
Before running this program, the C++ program must be executed with 10, 100, 1000, 10000 and 100000 timesteps.
This python file reads data from the different textfiles from the C++ program, and is implemented with the -03 compiler flags. 
Change the flags from True/False to run for the different files.

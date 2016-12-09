# Project5
This project solves the diffusion equation in one and two dimensions. By creating algorithms for the explicit and implicit Euler scheme, as well as for the implicit Crank-Nicolson scheme, the one-dimensional diffusion equation is solved by implementing these schemes in a C++ program. 
The two-dimensional diffusion equation is solved using the iterative Jacobi method and the iterative Gauss-Seidel method. The code is parallelized using OpenMP. This project also includes a program written in Python, which reads data from
output files before plotting the results. The programs are executed as follows:

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
This python file reads data from the different textfiles from the C++ program. Change the flags from True/False to run for the different files.

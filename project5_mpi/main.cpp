#define ARMA_NO_DEBUG
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <omp.h>
#include "solvers.h"

using namespace std;
using namespace arma;


int main(int argc, char * argv[]){

    // exit program if too few command line arguments
    if (argc < 3 ) {
        cout << "Error: bad usage! Please see the README.md file for command line arguments." << endl;
        exit(1);
    }

    // define variables
    int n = atoi( argv[1] );
    int timesteps = atoi( argv[2] );

    // Delta x variables to test stability -> change manually
    //double dx = 1.0 / (2.0 * n);
    //double dx = 1.0 / (4.0 * n);
    //double dx = 2.0 / n;
    //double dx = 3.0 / n;

    double dx = 1.0 / (n - 1);
    double dt = 0.25 * dx * dx;
    double alpha = dt / ( dx * dx );
    double tolerance = 1.0e-8;

    omp_set_num_threads(4);

    mat A(n ,n);
    mat A_init(n, n);

    // boundary conditions -> zeros at all endpoints
    for ( int i = 0.0; i < n; i ++ ){
        A(0, i) = 0.0;
        A(n-1, i) = 0.0;
        A(i, 0) = 0.0;
        A(i, n-1) = 0.0;
    }

    // initializing A_init
    for ( int i = 0; i < n; i++ ){
        for ( int j = 0; j < n; j++ ){
            A_init(i, j) = sin( 2 * M_PI * dx * i ) * sin( 2 * M_PI * dx * j );
        }
    }

    // start timer
    double wtime = omp_get_wtime ( );

    int iteration_counter;

    if ( atoi(argv[3]) == 1){

        // iterative Jacobi method

        for ( int t = 1; t < timesteps; t++ ){

            iteration_counter = JacobiSolver(n, tolerance, alpha, A, A_init);

            // set A_init equal to the calculated A, so that this becomes the previous step in next calculation
            for ( int i = 0; i < n; i ++ ){
                for ( int j = 0; j < n; j++ ){
                    A_init(i, j) = A(i, j);
                }
            }
        }
    }

    if ( atoi(argv[3]) == 2 ){

        // iterative Gauss-Seidel method

        for ( int t = 1; t < timesteps; t++ ){

            iteration_counter = GaussSeidelSolver( n, tolerance, alpha, A, A_init );

            for ( int i = 0; i < n; i ++ ){
                for ( int j = 0; j < n; j++ ){
                    A_init(i, j) = A(i, j);
                }
            }
        }
    }

    // exact solution for comparison
    double exact_solution;
    double sum = 0.0;
    for ( int t = 1; t <= timesteps; t++ ){
        for ( int i = 0; i < n; i++ ){
            for ( int j = 0; j < n; j++ ){
                exact_solution = sin( 2 * M_PI * dx * i ) * sin( 2 * M_PI * dx * j ) * exp( - 4 * M_PI * M_PI * dt * t );
                sum += fabs(A(i, j) - exact_solution);
            }
        }
    }

    // end timer
    wtime = omp_get_wtime ( ) - wtime;

    // write final solution matrix to file
    WriteToFile( n, A, argv[2] );

    cout << "Time used: " << wtime << endl;

    if ( atoi(argv[3]) == 1 ){
        cout << "Jacobi error is " << sum/( n * n ) << " in " << iteration_counter << " iterations" << endl;
    }

    if ( atoi(argv[3]) == 2 ){
        cout << "Gauss-Seidel error is " << sum/( n * n ) << " in " << iteration_counter << " iterations" << endl;
    }

    return 0;
}

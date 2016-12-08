#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <omp.h>
#include "solvers.h"

using namespace std;
using namespace arma;


int JacobiSolver( int n, double tolerance, double alpha, mat &A, mat &A_init ){

    int max_iterations = 100000;
    double D = 1 / (1 + 4 * alpha);   // matrix diagonal
    int i, j;

    mat A_old(n, n);


#pragma omp parallel for
        // fill A_old with alpha
        for ( int i = 1; i < n-1; i++ ){
            for ( int j = 1; j < n-1; j++ ){
                A_old(i, j) = alpha;
            }
        }

    // Jacobi iterative solver
    for ( int iterations = 0; iterations < max_iterations; iterations++ ){

#pragma omp parallel for
            for ( int i = 1; i < n-1; i++ ){
                for ( int j = 1; j < n-1; j++ ){
                    // implicit Euler scheme
                    A(i, j) = D * A_init(i, j) +
                            D * ( alpha * ( A_old(i+1, j) + A_old(i-1, j) + A_old(i, j+1) + A_old(i, j-1) ) );
                }
            }

        double sum = 0.0;
# pragma omp parallel default(shared) private(i, j) reduction(+:sum)
        {
# pragma omp for
            for ( int i = 0; i < n; i ++ ){
                for ( int j = 0; j < n; j++ ){
                    sum += ( A_old(i, j) - A(i, j) ) * ( A_old(i, j) - A(i, j) );
                    A_old(i, j) = A(i, j);
                }
            }
        }

        // convergence
        if ( sqrt(sum) < tolerance ){
            return iterations;
        }
    }

    cerr << "Maximum number of iterations reached without convergence" << endl;

    return max_iterations;
}

int GaussSeidelSolver( int n, double tolerance, double alpha, mat &A, mat &A_init ){

    int max_iterations = 100000;
    double D = 1 / ( 1 + 4 * alpha );   // matrix diagonal
    int i, j;

    mat A_old(n, n);


#pragma omp parallel for
        // fill A_old with alpha
        for ( int i = 1; i < n-1; i++ ){
            for ( int j = 1; j < n-1; j++ ){
                A_old(i, j) = alpha;
            }
        }

    // Gauss-Seidel iterative solver
    for ( int iterations = 0; iterations < max_iterations; iterations++ ) {  // k

#pragma omp parallel for
            for ( int i = 1; i < n-1; i++ ){
                for ( int j = 1; j < n-1; j++ ){
                    // implicit Euler scheme
                    A(i, j) = D * A_init(i, j) +
                            D * ( alpha * ( A_old(i+1, j) + A(i-1, j) + A_old(i, j+1) + A(i, j-1) ) );
                }
            }

        double sum = 0.0;
# pragma omp parallel default(shared) private(i, j) reduction(+:sum)
        {
# pragma omp for
            for ( int i = 0; i < n; i ++ ){
                for ( int j = 0; j < n; j++ ){
                    sum += ( A_old(i, j) - A(i, j) ) * ( A_old(i, j) - A(i, j) );
                    A_old(i, j) = A(i, j);
                }
            }
        }

        // convergence
        if ( sqrt(sum) < tolerance ){
            return iterations;
        }
    }

    cerr << "Maximum number of iterations reached without convergence" << endl;

    return max_iterations;
}

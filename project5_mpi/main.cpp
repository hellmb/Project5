#define ARMA_NO_DEBUG
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <omp.h>

using namespace std;
using namespace arma;


int JacobiSolver(int n, double tolerance, double alpha, mat &A, mat &A_init){

    int max_iterations = 100000;
    double D = 1 / (1 + 4 * alpha);   // matrix diagonal
    int i, j;

    mat A_old(n, n);


#pragma omp parallel for
        // fill A_old with alpha
        for (int i = 1; i < n-1; i++){
            for (int j = 1; j < n-1; j++){
                A_old(i, j) = alpha;
            }
        }

    // Jacobi iterative solver
    for (int iterations = 0; iterations < max_iterations; iterations ++){

#pragma omp parallel for
            for (int i = 1; i < n-1; i++){
                for (int j = 1; j < n-1; j++){
                    // implicit Euler scheme
                    A(i, j) = D * A_init(i, j) +
                            D * ( alpha * ( A_old(i+1, j) + A_old(i-1, j) + A_old(i, j+1) + A_old(i, j-1) ) );
                }
            }

        double sum = 0.0;
# pragma omp parallel default(shared) private(i, j) reduction(+:sum)
        {
# pragma omp for
            for (int i = 0; i < n; i ++){
                for (int j = 0; j < n; j++){
                    sum += ( A_old(i, j) - A(i, j) ) * ( A_old(i, j) - A(i, j) );
                    A_old(i, j) = A(i, j);
                }
            }
        }

        // convergence
        if (sqrt(sum) < tolerance){
            return iterations;
        }
    }

    cerr << "Maximum number of iterations reached without convergence" << endl;

    return max_iterations;
}



int GaussSeidelSolver(int n, double tolerance, double alpha, mat &A_jacobi, mat &A_init, mat &A){

    int max_iterations = 100000;
    double D = 1 / (1 + 4 * alpha);   // matrix diagonal
    int i, j;

    mat A_old(n, n);
    mat GS_prev(n, n);


#pragma omp parallel for
        // fill A_old with 1's
        for (int i = 1; i < n-1; i++){
            for (int j = 1; j < n-1; j++){
                A_old(i, j) = alpha;
            }
        }

    // do Jacobi iterative method once for first calculation of Gauss-Seidel
    for (int iterations = 0; iterations < max_iterations; iterations ++){

#pragma omp parallel for
            for (int i = 1; i < n-1; i++){
                for (int j = 1; j < n-1; j++){
                    // implicit Euler scheme
                    A_jacobi(i, j) = D * A_init(i, j) +
                            D * ( alpha * ( A_old(i+1, j) + A_old(i-1, j) + A_old(i, j+1) + A_old(i, j-1) ) );
                }
            }

        double sum = 0.0;
# pragma omp parallel default(shared) private(i, j) reduction(+:sum)
        {
# pragma omp for
            for (int i = 0; i < n; i ++){
                for (int j = 0; j < n; j++){
                    sum += ( A_old(i, j) - A_jacobi(i, j) ) * ( A_old(i, j) - A_jacobi(i, j) );
                    A_old(i, j) = A_jacobi(i, j);
                }
            }
        }

        // convergence
        if (sqrt(sum) < tolerance){
            break;
        }
    }

#pragma omp parallel for
    for (int i = 1; i < n-1; i++){
        for (int j = 1; j < n-1; j++){
            // set the previous value of of the GS solver to be the Jacobi matrix calculated above
            GS_prev(i, j) = A_jacobi(i, j);
        }
    }

    // Gauss-Seidel iterative method
    for (int iterations = 0; iterations < max_iterations; iterations++){

# pragma omp parallel for
        for (int i = 1; i < n-1; i++){
            for (int j = 1; j < n-1; j++){
                // implicit Euler scheme
//                A(i, j) = D * A_init(i, j) + D * GS_prev(i, j) +
//                        D * ( alpha * ( A_old(i+1, j) + A_old(i-1, j) + A_old(i, j+1) + A_old(i, j-1) ) );
                A(i, j) = D * A_init(i, j) + D * alpha * ( GS_prev(i-1, j) + GS_prev(i, j-1) + GS_prev(i+1, j) + GS_prev(i, j+1) ) +
                        D * alpha * ( A_old(i+1, j) + A_old(i-1, j) + A_old(i, j+1) + A_old(i, j-1)  );
            }
        }

        double sum;
# pragma omp parallel default(shared) private(i, j) reduction(+:sum)
        {
# pragma parallel for
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    sum += ( A_old(i, j) - A(i, j) ) * ( A_old(i, j) - A(i, j) );
                    //sum += ( A(i, j) - A_old(i, j) ) * ( A(i, j) - A_old(i, j) );
                    //A_old(i, j) = GS_prev(i, j);
                    //GS_prev(i, j) = A(i, j);
                    A_old(i, j) = A(i, j);
                }
            }
        }

        // convergence
        if (sqrt(sum) < tolerance){
            return iterations;
        }
    }

    cerr << "Maximum number of iterations reached without convergence" << endl;

    return max_iterations;
}



int main(int argc, char * argv[]){

    // exit program if too few command line arguments
    if (argc < 3 ) {
        cout << "Error: bad usage! Please set the matrix dimension and preferred number of timesteps." << endl;
        exit(1);
    }

    // define variables
    int n = atoi(argv[1]);
    int timesteps = atoi(argv[2]);
    double dx = 1.0/(n*n);
    double dt = 0.25 * dx * dx;
    double alpha = dt / (dx * dx);
    double tolerance = 1.0e-14;

    omp_set_num_threads(4);

    mat A(n ,n);
    mat A_init(n, n);
    mat A_jacobi(n, n);

    // boundary conditions -> zeros at all endpoints
    for (int i = 0.0; i < n; i ++){
        A(0, i) = 0.0;
        A(n-1, i) = 0.0;
        A(i, 0) = 0.0;
        A(i, n-1) = 0.0;
    }

    // initializing A_jacobi
    for (int i = 0; i < n; i ++){
        for (int j = 0; j < n; j++){
            A_jacobi(i, j) = A(i, j);
        }
    }


    double wtime = omp_get_wtime ( );

    int iteration_counter;

    if (atoi(argv[3]) == 1){

        for (int t = 1; t < timesteps; t++){

            iteration_counter = JacobiSolver(n, tolerance, alpha, A, A_init);

            // set A_init equal to the calculated A, so that this becomes the previous step in next calculation
            for (int i = 0; i < n; i ++){
                for (int j = 0; j < n; j++){
                    A_init(i, j) = A(i, j);
                }
            }
        }
    }

    if (atoi(argv[3]) == 2){

        for (int t = 1; t < timesteps; t++){

            iteration_counter = GaussSeidelSolver(n, tolerance, alpha, A_jacobi, A_init, A);

            for (int i = 0; i < n; i ++){
                for (int j = 0; j < n; j++){
                    A_init(i, j) = A(i, j);
                }
            }
        }
    }

    //A_jacobi.print("A jacobi: ");

    //A.print("A: ");


    // exact solution for comparison
    double exact_solution;
    double sum = 0.0;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            exact_solution = -sin(M_PI * dx * i) * sin(M_PI * dx * j);
            sum += fabs(A(i, j) - exact_solution);
        }
    }

    wtime = omp_get_wtime ( ) - wtime;

    cout << "Time used: " << wtime << endl;
    cout << "Jacobi error is " << sum/(n*n) << " in " << iteration_counter << " iterations" << endl;

    return 0;
}


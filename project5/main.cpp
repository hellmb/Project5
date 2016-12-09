#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <solvers.h>

using namespace std;
using namespace arma;


int main ( int argc, char* argv[] ){

    // exit program if too few command line arguments
    if (argc < 2 ) {
        cout << "Error: bad usage! Please see the README.md file for command line arguments.";
        exit(1);
    }

    // define variables
    int n = 100;
    int timesteps = atoi(argv[1]);

    double delta_x = 0.1;
    double delta_t = ( delta_x * delta_x ) / 2.0;

    // different alpha value for Euler and Crank-Nicolson schemes
    double alpha = delta_t / (delta_x * delta_x);
    double alpha_CN = delta_t / (2.0 * delta_x * delta_x);


    vec V_old(n);
    vec V_new(n);

    // boundary conditions
    V_old(0) = 0.0;     V_old(n-1) = 1.0;
    V_new(0) = 0.0;     V_new(n-1) = 1.0;

    if ( atoi(argv[2]) == 1 ){

        // exlicit forward Euler scheme
        explicitForwardEuler( n, timesteps, V_old, V_new, alpha, argv[1] );

    }

    if ( atoi(argv[2]) == 2 ) {

        // implicit backward Euler scheme

        double a = - alpha;
        double c = - alpha;
        double b = 1.0 + 2.0 * alpha;

        for ( int t = 1; t <= timesteps; t++ ){

            TridiagonalSolver( a, b, c, n, V_new, V_old );

            V_old = V_new;
        }
    }

    if ( atoi(argv[2]) == 3 ){

        // implicit Crank-Nicolson scheme

        double a = - alpha_CN;
        double c = - alpha_CN;
        double b = 2.0 + 2.0 * alpha_CN;

        for ( int t = 1; t <= timesteps; t++ ){

            for ( int i = 1; i < n-1; i++ ){
                // explicit scheme
                V_old(i) = alpha_CN * V_new(i-1) + (2 - 2 * alpha_CN) * V_new(i) + alpha_CN * V_new(i+1);

            }

            TridiagonalSolver( a, b, c, n, V_new, V_old );

        }

        // write file for different timsteps
        string implicit_CN ("../files_project5/implicit_CN");
        implicit_CN += argv[1];
        implicit_CN += ".txt";

        ofstream myfile3;
        myfile3.open( implicit_CN );

        for ( int i = 0; i < n; i++ ){
            myfile3 << V_new(i) << endl;
        }
    }

    if ( atoi(argv[2]) == 4 ){

        // explicit Euler scheme for unstable solutions
        explicitForwardEulerUnstable( n, timesteps, V_old, V_new, argv[1] );
    }


    return 0;
}





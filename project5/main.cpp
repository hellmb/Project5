/*
 Input arguments:
 - argv[1]: 0.1 or 0.01 -> delta_x value
 - argv[2]: number of timesteps
 - argv[3]: 1, 2 or 3   -> deciding whcich scheme to solve
*/
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;


// tridiagnal solver function
void TridiagonalSolver(double a, double b, double c, int n, vec &V_new, vec &V_old){

    // define constants and vectors
    double constant = a * c;
    vec k = zeros<vec>(n);
    vec diag = zeros<vec>(n);

    // initialize diagonal vector
    diag(0) = b;

    // forward substitution
    for (int i = 1; i < n; i++){
        diag(i) = b - constant / diag(i-1);
        k(i) = a / diag(i);
        V_old(i) = V_old(i) - V_old(i-1) * k(i-1);
    }

    // backward substitution
    for (int i = n-2; i > 0; i--){
        V_new(i) = ( V_old(i) - c * V_new(i+1) ) / diag(i);
    }
}

void explicitForwardEuler(int n, int timesteps, vec &V_old, vec &V_new, double alpha, string timestepsString) {

    for (int t = 1; t <= timesteps; t++){
        for (int i = 1; i < n-1; i++){
            V_new(i) = alpha * V_old(i-1) + (1.0 - 2.0 * alpha) * V_old(i) + alpha * V_old(i+1);
        }
        V_old = V_new;
    }
    // write file for different timsteps
    string explicit_euler ("../files_project5/explicit_euler");
    explicit_euler += timestepsString;
    explicit_euler += ".txt";

    ofstream myfile1;
    myfile1.open(explicit_euler);

    for (int i = 0; i < n; i++){
        myfile1 << V_new(i) << endl;
    }
}

int main(int argc, char* argv[]){

    // exit program if too few command line arguments
    if (argc <= 3 ) {
        cout << "Error: bad usage! Please see the README.mp file for command line arguments.";
        exit(1);
    }

    // define variables

    int n = 100;
    int timesteps = atoi(argv[2]);

    double delta_x = atof(argv[1]);
    double delta_t = ( delta_x * delta_x ) / 2.0;

    // test unstable solution for the explicit scheme
    //double delta_t = 1.0;

    double alpha = delta_t / (delta_x * delta_x);


    vec V_old = zeros<vec>(n);
    vec V_new = zeros<vec>(n);

    // boundary conditions
    V_old(0) = 0.0;     V_old(n-1) = 1.0;
    V_new(0) = 0.0;     V_new(n-1) = 1.0;

    if (atoi(argv[3]) == 1){
        // exlicit forward Euler scheme

        explicitForwardEuler(n, timesteps, V_old, V_new, alpha, argv[2]);

    }

    if (atoi(argv[3]) == 2){
        // implicit backward Euler scheme

        double a = - alpha;
        double c = - alpha;
        double b = 1.0 + 2.0 * alpha;

        for (int t = 1; t <= timesteps; t++){

            TridiagonalSolver(a, b, c, n, V_new, V_old);

            V_old = V_new;
        }

        V_new.print("V_new: ");

        // write file for different timsteps
        string implicit_euler ("../files_project5/implicit_euler");
        implicit_euler += argv[2];
        implicit_euler += ".txt";

        ofstream myfile3;
        myfile3.open(implicit_euler);

        for (int i = 0; i < n; i++){
            myfile3 << V_new(i) << endl;
        }

    }

    if (atoi(argv[3]) == 3){
        // implicit Crank-Nicolson scheme

        double a = - alpha;
        double c = - alpha;
        double b = 2.0 + 2.0 * alpha;

        vec V_r = zeros<vec>(n);
        V_r(0) = 0.0;
        V_r(n-1) = 1.0;

        vec V_CN = zeros<vec>(n);
        V_CN(0) = 0.0;
        V_CN(n-1) = 1.0;

        for (int t = 1; t <= timesteps; t++){

            for (int i = 1; i < n-1; i++){
                // explicit scheme
                V_r(i) = alpha * V_CN(i-1) + (2 - 2 * alpha) * V_CN(i) + alpha * V_CN(i+1);

            }

            TridiagonalSolver(a, b, c, n, V_CN, V_r);

        }

        // write file for different timsteps
        string implicit_CN ("../files_project5/implicit_CN");
        implicit_CN += argv[2];
        implicit_CN += ".txt";

        ofstream myfile3;
        myfile3.open(implicit_CN);

        for (int i = 0; i < n; i++){
            myfile3 << V_CN(i) << endl;
        }
    }

    if (atoi(argv[3]) == 4){
        // Jacobi iterative method for explicit scheme

        mat b = zeros<mat>(n, n);
        mat b_temp = zeros<mat>(n, n);
        vec x = zeros<vec>(n);

        for (int i = 0; i < n; i++){
            x(i) = delta_x * delta_x * i;
        }

        // boundary
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (j == n-1){
                    b(i, j) = 100.0 * sin( (M_PI * x(i))/x(n-1) );
                }
                if (j == 0){
                    b(i, j) = -100.0 * sin( (M_PI * x(i))/x(n-1) );
                }
            }
        }
         // solve Jacobi iterative method over timesteps

        // loop over time
        for (int t = 1; t < timesteps; t++){

            int iterations = 0;
            int max_iterations = 100;
            double tolerance = 0.00001;
            double difference = tolerance + 1;

            while (iterations <= max_iterations && difference > tolerance){

                b_temp = b;
                difference = 0;
                for (int i = 1; i < n-1; i++){
                    for (int j = 1; j < n-1; j++){

                        b(i, j) = 0.25 * ( b_temp(i, j+1) + b_temp(i, j-1) + b_temp(i+1, j) + b_temp(i-1, j) );
                        // implicit Euler scheme
                        //b(i, j) = (1.0 + 4.0 * alpha) * b_temp(i, j) - alpha * (b_temp(i+1, j) + b_temp(i-1, j) + b_temp(i, j+1) + b_temp(i, j-1));
                        difference += fabs( sum(b(i, j) - b_temp(i, j)) );
                    }
                }
                difference /= n * n;
                iterations += 1;
                //cout << "Difference: " << difference << endl;
            }
        }
        b.print();




//        // defining two-dimensional matrix
//        mat A = zeros<mat>(n, n);
//        mat D = zeros<mat>(n, n);
//        mat L = zeros<mat>(n, n);
//        mat u = zeros<mat>(n, n);

//        vec temp(n);
//        for (int i = 1; i < n-1; i++){
//            temp(i) = -1.0;
//        }

//        // fill matrices
//        for (int i = 0; i < n; i++){
//            for (int j = 0; j < n; j++){

//                if (i == 0){
//                    u(i, j) = temp(j);
//                }
//                if (j == n-1){
//                    u(i, j) = temp(i);
//                }

//                if (j == 0){
//                    L(i, j) = temp(i);
//                }
//                if (i == n-1){
//                    L(i, j) = temp(j);
//                }

//                if (i == j){
//                    D(i, j) = 4.0;
//                }
//            }
//        }

//        A = D + L + u;

//        //A.print("A: ");

//        vec b(n);
//        for (int i = 0; i < n; i++){
//            // this is wrong
//            b(i) = delta_x * delta_x;
//        }

//        //b.print("b: ");

//        vec sol_new = zeros<vec>(n);
//        // initial solution guessed at 0
//        vec sol_old = randu<vec>(n);


//        // create a loop over time here?
//        //for (int t = 1; t < timesteps; t++){

//        // start iteration algorithm
//        int iterations = 0;
//        int max_iterations = 1000;
//        double tolerance = 0.00001;
//        double difference = tolerance + 1.0;

//            while ( iterations <= max_iterations && difference > tolerance ) {

//                difference = 0;

//                sol_new = 0.25 * (b - (L+u) * sol_old);
//                iterations += 1;
//                difference += fabs(sum(sol_new - sol_old)/n);

//                sol_old = sol_new;
//                //cout << difference << endl;

//            }
//            //difference /= pow(n, 2.0);
//            cout << difference << endl;
//        //}

//        vec solution = zeros<vec>(n);
//        solution = sol_old;

//        solution.print("Solution: ");
   }


    return 0;
}


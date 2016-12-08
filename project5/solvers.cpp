#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <solvers.h>

using namespace std;
using namespace arma;

// explicit Euler scheme
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

// explicit Euler scheme for unstable solutions
void explicitForwardEulerUnstable(int n, int timesteps, vec &V_old, vec &V_new, string timestepsString) {

    double alpha_unstable = 1.0;

    for (int t = 1; t <= timesteps; t++){
        for (int i = 1; i < n-1; i++){
            V_new(i) = alpha_unstable * V_old(i-1) + (1.0 - 2.0 * alpha_unstable) * V_old(i) + alpha_unstable * V_old(i+1);
        }
        V_old = V_new;
    }
    // write file for different timsteps
    string explicit_euler ("../files_project5/explicit_euler_unstable");
    explicit_euler += timestepsString;
    explicit_euler += ".txt";

    ofstream myfile4;
    myfile4.open(explicit_euler);

    for (int i = 0; i < n; i++){
        myfile4 << V_new(i) << endl;
    }
}

// tridiagnal solver for the implcit schemes
void TridiagonalSolver(double a, double b, double c, int n, vec &V_new, vec &V_old){

    // define constants and vectors
    double constant = a * c;
    vec k(n);
    vec diag(n);

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

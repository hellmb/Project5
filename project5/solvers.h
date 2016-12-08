#ifndef SOLVERS_H
#define SOLVERS_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void explicitForwardEuler(int n, int timesteps, vec &V_old, vec &V_new, double alpha, string timestepsString) ;
void explicitForwardEulerUnstable(int n, int timesteps, vec &V_old, vec &V_new, string timestepsString)  ;
void TridiagonalSolver(double a, double b, double c, int n, vec &V_new, vec &V_old) ;

#endif // SOLVERS_H

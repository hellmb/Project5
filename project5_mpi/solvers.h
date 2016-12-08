#ifndef SOLVERS_H
#define SOLVERS_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int JacobiSolver( int n, double tolerance, double alpha, mat &A, mat &A_init ) ;
int GaussSeidelSolver( int n, double tolerance, double alpha, mat &A, mat &A_init ) ;
void WriteToFile( int n, mat &A, string timestepString ) ;

#endif // SOLVERS_H

#pragma once
#include <vector>
#include <cmath>

namespace MATH60082
{
/* 
ON INPUT:
a, b and c -- are the diagonals of the matrix
rhs        -- is the right hand side
x          -- is the initial guess
iterMax    -- is maximum iterations
tol        -- is the tolerance level
omega      -- is the relaxation parameter 
sor        -- not used
ON OUTPUT:
a, b, c, rhs        -- unchanged
x                   -- solution to Ax=b
iterMax, tol, omega -- unchanged
sor                 -- number of iterations to converge
*/
void sorSolve(const std::vector<double> &a,const std::vector<double> &b,const std::vector<double> &c,const std::vector<double> &rhs,
std::vector<double> &x,int iterMax,double tol,double omega,int &sor );

}

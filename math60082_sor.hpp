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
std::vector<double> &x,int iterMax,double tol,double omega,int &sor )
{
    // assumes vectors a,b,c,d,rhs and x are same size (doesn't check)
    int n=a.size()-1;
    // sor loop
  for(sor=0;sor<iterMax;sor++)
  {
    double error=0.;
    // implement sor in here
    {
      double y = (rhs[0] - c[0]*x[1])/b[0];
      x[0] = x[0] + omega*(y-x[0]); 
    }
    for(int j=1;j<n;j++)
    {
      double y = (rhs[j] - a[j]*x[j-1] - c[j]*x[j+1])/b[j];
      x[j] = x[j] + omega*(y-x[j]); 
    }
    {
      double y = (rhs[n] - a[n]*x[n-1])/b[n];
      x[n] = x[n] + omega*(y-x[n]); 
    }
    // calculate residual norm ||r|| as sum of absolute values
    error += std::fabs(rhs[0] - b[0]*x[0] - c[0]*x[1]);
    for(int j=1;j<n;j++) 
        error += std::fabs(rhs[j] - a[j]*x[j-1] - b[j]*x[j] - c[j]*x[j+1]);
    error += std::fabs(rhs[n] - a[n]*x[n-1] - b[n]*x[n]);
    // make an exit condition when solution found
    if(error<tol)
      break;
  }
}
}

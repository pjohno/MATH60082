#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

namespace MATH60082
{
/* 
ON INPUT:
S0          -- initial stock price
X           -- exercise (strike) price
T           -- Time to expiry (years)
r           -- interest rate (per annum)
sigma       -- volatility (per annum^12)
iMax        -- number of time steps
jMax        -- number of space steps
SMax        -- Maximum value of S
omega       -- is the relaxation parameter 
tol         -- is the tolerance level
iterMax     -- is maximum iterations
ON OUTPUT:
return      -- the value of a European put option at S=S0, t=0
*/
double europeanPut_CN(double S0,double X,double T,double r,double sigma,int iMax,int jMax,double SMax,double omega,double tol,int iterMax);

double normalDistribution(double x);

// return the value of a put option using the black scholes formula
double europeanPut_exact(double S,double X,double T,double r,double sigma);

}

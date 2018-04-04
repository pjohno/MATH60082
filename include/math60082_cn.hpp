#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

namespace MATH60082
{
/** 
@brief Return the value of a European Put option using the Crank-Nicolson method with SOR solver
ON INPUT:
@param S0          -- initial stock price
@param X           -- exercise (strike) price
@param T           -- Time to expiry (years)
@param r           -- interest rate (per annum)
@param sigma       -- volatility (per annum^12)
@param iMax        -- number of time steps
@param jMax        -- number of space steps
@param SMax        -- Maximum value of S
@param omega       -- is the relaxation parameter 
@param tol         -- is the tolerance level
@param iterMax     -- is maximum iterations
ON OUTPUT:
@return the value of a European put option at \f$ S=S0 \f$, \f$ t=0 \f$
**/
double europeanPut_CN(double S0,double X,double T,double r,double sigma,int iMax,int jMax,double SMax,double omega,double tol,int iterMax);

/** 
@brief Return the value of an American Put option using the Crank-Nicolson method with PSOR solver
ON INPUT:
@param S0          -- initial stock price
@param X           -- exercise (strike) price
@param T           -- Time to expiry (years)
@param r           -- interest rate (per annum)
@param sigma       -- volatility (per annum^12)
@param iMax        -- number of time steps
@param jMax        -- number of space steps
@param SMax        -- Maximum value of S
@param omega       -- is the relaxation parameter 
@param tol         -- is the tolerance level
@param iterMax     -- is maximum iterations
ON OUTPUT:
@return the value of a American put option at \f$ S=S0 \f$, \f$ t=0 \f$
**/
double americanPut_CN(double S0,double X,double T,double r,double sigma,int iMax,int jMax,double SMax,double omega,double tol,int iterMax);

/* 
ON INPUT:
S0          -- initial stock price
*/
double normalDistribution(double x);

// return the value of a put option using the black scholes formula
double europeanPut_exact(double S,double X,double T,double r,double sigma);

}

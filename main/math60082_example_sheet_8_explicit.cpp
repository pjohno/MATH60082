#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

/* Template code for the Explicit Finite Difference
*  American put option with
*     V(S,T) = max( X-S , 0 )
*  and
*      dV/dt + 1/2 sigma^2 S^2 d^2V/dS^2 + r S dV/dS - r V = 0
*  with American condition
*         V >= X - S
*  At S=0 we have
*     V = X
*  As S->infty, assume V->0.
*/
double AmericanPut_explicit(double S0, double X, double T, double r, double sigma, int iMax, int jMax, double S_max)
{
	// declare and initialise local variables (ds,dt)
	double dS = S_max / jMax;
	double dt = T / iMax;
	// create storage for the stock price and option price (old and new)
	vector<double> S(jMax + 1), vOld(jMax + 1), vNew(jMax + 1);
	// setup and initialise the stock price 
	for (int j = 0; j <= jMax; j++)
	{
		S[j] = j*dS;
	}
	// setup and initialise the final conditions on the option price 
	for (int j = 0; j <= jMax; j++)
	{
		vOld[j] = max(X - S[j],0.);
		vNew[j] = max(X - S[j],0.);
	}
	// loop through time levels, setting the option price at each grid point, and also on the boundaries
	for (int i = iMax - 1; i >= 0; i--)
	{
		// apply boundary condition at S=0
		vNew[0] = X;
		for (int j = 1; j <= jMax - 1; j++)
		{
			double A, B, C;
			// matches PDE at the top
			A = 0.5*sigma*sigma*j*j*dt + 0.5*r*j*dt;
			B = 1. - sigma*sigma*j*j*dt;
			C = 0.5*sigma*sigma*j*j*dt - 0.5*r*j*dt;
                        // European or continuation value
			vNew[j] = 1. / (1. + r*dt)*(A*vOld[j + 1] + B*vOld[j] + C*vOld[j - 1]);
			// American condition says
			// V >= X - S 
			// so take 
			// max( X - S , V )
			vNew[j] = max( X - S[j] , vNew[j] );
		}
		// apply boundary condition at S=S_max
		vNew[jMax] = 0.;
		// set old values to new
		vOld = vNew;
	}

	// get j* such that S_0 \in [ j*dS , (j*+1)dS ]
	int jstar;
	jstar = S0 / dS;
	double sum = 0.;
	// run 2 point Lagrange polynomial interpolation
	sum = sum + (S0 - S[jstar + 1]) / (S[jstar] - S[jstar + 1])*vNew[jstar];
	sum = sum + (S0 - S[jstar]) / (S[jstar + 1] - S[jstar])*vNew[jstar + 1];
	return sum;
}

int main()
{
	// declare and initialise Black Scholes parameters
	double S0 = 9.735, X = 10., T = 1., r = 0.05, sigma = 0.4;
	// declare and initialise grid paramaters 
	int iMax = 400, jMax = 40;
	cout << AmericanPut_explicit(S0, X, T, r, sigma, iMax, jMax, 5.*X) << endl;
	/*
         * OUTPUT >>
         * 1.47058
	*/
}

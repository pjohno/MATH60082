#include "math60082_monteCarloLeastSquares.hpp"
#include <iostream>
using namespace std;

int main()
{
    // number of exercise dates
    int timeSteps=10; 
    // some parameters in the problem
    double maturity=1.,interestRate=0.05,strikePrice=2.,S0=1.973,sigma=0.4;
    
    // number of paths n
    int n=100;
    // degree of fitting polynomial
    int M_degree=4;
    
    for(int i=1;i<=1000;i++)
    {
        int N=i*n;
        cout << N << " " << MATH60082::monteCarloLeastSquares(S0,maturity,strikePrice,interestRate,sigma,timeSteps,M_degree,N) << endl;
    }
}

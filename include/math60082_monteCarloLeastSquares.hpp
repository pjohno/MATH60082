#ifndef _COMPFIN_MONTECARLOLEASTSQUARES_HPP_
#define _COMPFIN_MONTECARLOLEASTSQUARES_HPP_

namespace MATH60082
{
    
    double monteCarloLeastSquares(double S0,double maturity,double strikePrice,
                                  double interestRate,double sigma,int timeSteps,int M_degree,int n);
}

#endif 

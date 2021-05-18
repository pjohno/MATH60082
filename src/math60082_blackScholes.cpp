#include "math60082_blackScholes.hpp"
#include <boost/math/distributions.hpp>
#include <cmath>
 
namespace MATH60082{
    
    double blackScholesPutOption(double S,double X,double T,double r,double d,double sigma)
    { 
        static boost::math::normal_distribution<> nd;
        if(S<1.e-14)return X*exp(-r*T);
        if(T<1.e-14){if(S>X)return 0.;else return X-S;}
        double d1=(log(S/X)+(r-d+sigma*sigma/2.)*T)/sigma/sqrt(T);
        double d2=d1-sigma*sqrt(T);
        return X*exp(-r*T)*cdf(nd,-d2) - S*exp(-d*T)*cdf(nd,-d1);
    }
    
    double blackScholesCallOption(double S,double X,double T,double r,double d,double sigma)
    { 
        static boost::math::normal_distribution<> nd;
        if(S<1.e-14)return 0.;
        if(T<1.e-14){if(S<X)return 0.;else return S-X;}
        double d1=(log(S/X)+(r-d+sigma*sigma/2.)*T)/sigma/sqrt(T);
        double d2=d1-sigma*sqrt(T);
        return S*exp(-d*T)*cdf(nd,d1) - X*exp(-r*T)*cdf(nd,d2);
    }
    
}

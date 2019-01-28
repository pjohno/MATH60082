#ifndef _COMPFIN_LEASTSQUARESFIT_HPP_
#define _COMPFIN_LEASTSQUARESFIT_HPP_
#include "math60082_cheb.hpp"
#include "math60082_data.hpp"
// standard libraries
#include <algorithm>

namespace MATH60082
{
    
    
    class LeastSquaresFit
    {
        
        // we will generate a Chebyshev polynomial
        // of degree k
        double chiSq;
        ChebyshevPolynomial P;
        std::vector<std::vector<double>> C;
        double xTotal,xRange;
        
    public:
        
        void generateFit(int k,std::vector<DataPoint> &data);
        
        // return the fitted value of degree n
        double operator()(double x,unsigned int n) const 
        {
            unsigned int i=std::min(n,P.size()-1);
            return P(x,i);
        }
        // return value using highest degree
        double operator()(double x) const 
        {
            return P(x);
        }
        
        unsigned int size() const { return P.size(); }
        
        friend std::ostream& operator<<(std::ostream& output,const LeastSquaresFit& fit);
        
    };
    
}

#endif

#ifndef _COMPFIN_CHEBYSHEVPOLYNOMIAL_HPP_
#define _COMPFIN_CHEBYSHEVPOLYNOMIAL_HPP_

// c++ stuff
#include <iostream>
#include <cmath>
#include <vector>

//#define _PVJ_DEBUG_

namespace MATH60082 
{
    // 
    class ChebyshevPolynomial
    {
        // store the polynomial coefficients
        std::vector<double> cs; 
        double a,b;
        
    public:
        // default constructor
        ChebyshevPolynomial(){}
        // or vector
        ChebyshevPolynomial(const std::vector<double> &cs,double a,double b):cs(cs),a(a),b(b){}
        // or another ChebyshevPolynomial
        ChebyshevPolynomial(const ChebyshevPolynomial &C):cs(C.cs),a(C.a),b(C.b){}
        
        // equality can be dealt with like this for ChebyshevPolynomial
        ChebyshevPolynomial& operator=(const ChebyshevPolynomial &C)
        {
            cs=C.cs;
            a=C.a;
            b=C.b;
            return *this;
        }
        void init(int p,double a,double b);
        // return the fitted value of degree n
        double& operator[](int i)
        {
            return cs[i];
        }
        // return the fitted value of degree n
        double operator[](int i) const 
        {
            return cs[i];
        }
        // return the value of the chebyshev poly with error check
        double operator()(double x) const 
        {
            
            size_t i;
            double d1 = 0.0;
            double d2 = 0.0;
            
            double y = (2.0 * x - a - b) / (b - a);
            double y2 = 2.0 * y;
            
            for (i = cs.size()-1; i >= 1; i--)
            {
                double temp = d1;
                d1 = y2 * d1 - d2 + cs[i];
                d2 = temp;
            }
            
            return y * d1 - d2 +  cs[0];
        }
        // return the value of the chebyshev poly with error check
        double operator()(double x,int k) const 
        {
            
            size_t i;
            double d1 = 0.0;
            double d2 = 0.0;
            
            double y = (2.0 * x - a - b) / (b - a);
            double y2 = 2.0 * y;
            
            for (i = k; i >= 1; i--)
            {
                double temp = d1;
                d1 = y2 * d1 - d2 + cs[i];
                d2 = temp;
            }
            
            return y * d1 - d2 + cs[0];
        }
        
        unsigned int size() const {return cs.size();}
        
        // print the polynomial
        friend std::ostream& operator<<(std::ostream& output,const ChebyshevPolynomial& C);
        
    };
    
}

#endif

#include "math60082_cheb.hpp"
using namespace std;

namespace MATH60082{
    
    void ChebyshevPolynomial::init(int p,double a_,double b_)
    {
        cs = std::vector<double>(p);
        a = a_;
        b = b_;
    }
    
    std::ostream& operator<<(std::ostream& output,const ChebyshevPolynomial& C)
    {
        output << " Chebyshev polynomial degree ";
        output << C.cs.size()-1 << " :: { ";
        for(auto csi : C.cs)
            output << csi << " ";
        output << "}";
        return output;
    }
    
}

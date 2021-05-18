#include "math60082_stencils.hpp"
#include <algorithm>

namespace MATH60082
{
    
    std::ostream& operator<<(std::ostream& output,const Stencil& X)
    {
        output << " | a " << X.a << " | b " << X.b << " | c " << X.c <<" | ";
        return output;
    }
    std::ostream& operator<<(std::ostream& output,const FullStencil& X)
    {
        output << " | a " << X.a << " | b " << X.b << " | c " << X.c <<" | d "<< X.d <<" | ";
        return output;
    }
    
    std::ostream& operator<<(std::ostream& output,const VectorStencil& X)
    {
        for(unsigned int i=0;i<X.size();i++)
        {
            output << " i=" << i << " " << X[i] <<"\n";
        }
        return output;
    }
    
    int VectorStencil::createStencils(const std::function<double(double,double)> &alpha,const std::function<double(double,double)> &beta,const std::vector<double>& x,const double &t)
    {
        
        S.resize(x.size());
        double a_coeff,b_coeff;
        // default condition at x_min
        {
            size_t i=0;
            a_coeff = std::max(0.,alpha(x[i],t));
            S[i].a = 0.;
            S[i].b = -a_coeff/(x[i+1]-x[i]);
            S[i].c = a_coeff/(x[i+1]-x[i]);
        }
        // midpoints
        for(size_t i=1;i<x.size()-1;i++)
        {
            a_coeff = alpha(x[i],t);
            b_coeff = beta(x[i],t);
            b_coeff = 0.5*b_coeff*b_coeff/(x[i+1]-x[i])/(x[i]-x[i-1]);
            S[i].a = b_coeff;
            S[i].b = -2*b_coeff;
            S[i].c = b_coeff;
            if(a_coeff/(x[i+1]-x[i-1])>b_coeff)
            {
                S[i].b = S[i].b-a_coeff/(x[i+1]-x[i]);
                S[i].c = S[i].c+a_coeff/(x[i+1]-x[i]);
            }
            else if(a_coeff/(x[i+1]-x[i-1])<-b_coeff)
            {
                S[i].a = S[i].a-a_coeff/(x[i]-x[i-1]);
                S[i].b = S[i].b+a_coeff/(x[i]-x[i-1]);
            }
            else
            {
                S[i].a = S[i].a-a_coeff/(x[i+1]-x[i-1]);
                S[i].c = S[i].c+a_coeff/(x[i+1]-x[i-1]);
            }
            
        }
        
        // default condition at x_max
        {
            size_t i=x.size()-1;
            a_coeff = std::min(0.,alpha(x[i],t));
            S[i].a = -a_coeff/(x[i]-x[i-1]);
            S[i].b = a_coeff/(x[i]-x[i-1]);
            S[i].c = 0.;
            
        }
        
        return 1;
        
    }
    
    
}

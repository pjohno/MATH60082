#include "math60082_tridag.hpp"



namespace MATH60082
{
    /* 
     *    ON INPUT:
     *    a, b and c -- are the diagonals of the matrix
     *    rhs        -- is the right hand side
     *    ON OUTPUT:
     *    a, b, c             -- unchanged
     *    rhs                 -- solution to Ax=b
     */
    void thomasSolve(const std::vector<double> &a,const std::vector<double> &b_,const std::vector<double> &c,std::vector<double> &rhs)
    {
        int n=a.size();
        std::vector<double> b(n);
        // initial first value of b
        b[0]=b_[0];
        for(int j=1;j<n;j++)
        {
            b[j]=b_[j]-c[j-1]*a[j]/b[j-1];  
            rhs[j]=rhs[j]-rhs[j-1]*a[j]/b[j-1];
        }
        // calculate solution
        rhs[n-1]=rhs[n-1]/b[n-1];
        for(int j=n-2;j>=0;j--)
            rhs[j]=(rhs[j]-c[j]*rhs[j+1])/b[j];
    }
    
}

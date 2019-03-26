#pragma once
#include <vector>
#include <cmath>

namespace MATH60082
{
    /* 
     * ON INPUT:
     * a, b and c -- are the diagonals of the matrix
     * rhs        -- is the right hand side
     * ON OUTPUT:
     * a, b, c             -- unchanged
     * rhs                 -- solution to Ax=b
     */
    void thomasSolve(const std::vector<double> &a,const std::vector<double> &b,const std::vector<double> &c,std::vector<double> &rhs);
    
}


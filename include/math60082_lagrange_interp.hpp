#pragma once
#include <vector>

namespace MATH60082
{
    // A generic lagrange interpolation function
    double lagrangeInterpolation(const std::vector<double>& y,const std::vector<double>& x,double x0,unsigned int n=4);
}

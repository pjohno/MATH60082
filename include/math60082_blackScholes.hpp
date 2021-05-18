#pragma once
#include <iostream>
#include <vector>
#include <functional>

namespace MATH60082{

    
    double blackScholesPutOption(double S,double X,double T,double r,double d,double sigma);
    double blackScholesCallOption(double S,double X,double T,double r,double d,double sigma);
    
}

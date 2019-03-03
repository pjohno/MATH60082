#ifndef _COMPFIN_RANDOMPATHS_HPP_
#define _COMPFIN_RANDOMPATHS_HPP_

// c++ stuff
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

namespace MATH60082 
{
    std::vector< std::vector<double> > generateRandomPaths(int timesteps,double S0,double dt,double mu,double sigma,int N);
    
    std::vector< std::vector<double> > generateLongstaffPaths(int &timesteps,double &dt,int &N);
}

#endif 

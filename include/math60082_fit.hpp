#pragma once
#include "math60082_data.hpp"
#include <vector>

namespace MATH60082
{
    int runFitGSL(std::vector<DataPoint> &data,int p);
    int runFitEigen(std::vector<DataPoint> &data,int p);
}

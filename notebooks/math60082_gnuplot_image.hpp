#include <vector>
#include <iostream>
#include <iosfwd>
#include <sstream>

#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/process.hpp>
#include <boost/chrono.hpp>
#include "xeus/xjson.hpp"

namespace MATH60082{
    
    
    // allow for output from gnuplot to be plotted
    xeus::xjson mime_bundle_repr(const gnuplotImage& i)
    {
        auto bundle = xeus::xjson::object();
        bundle["image/png"] = encode64(i.imageText);
        return bundle;
    }
    
    // gnuplot functions
    GnuplotWidget G;
    gnuplotImage figure;
    // storage for plots
    std::vector<double> x;
    std::vector<double> y;
    std::vector<std::vector<double>> z;
}


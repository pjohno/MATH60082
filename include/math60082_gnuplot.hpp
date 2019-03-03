#pragma once
#include <sstream>
#include <vector>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/algorithm/string.hpp>

namespace MATH60082
{
  std::string decode64(const std::string &val) {
    using namespace boost::archive::iterators;
    using It = transform_width<binary_from_base64<std::string::const_iterator>, 8, 6>;
    return boost::algorithm::trim_right_copy_if(std::string(It(std::begin(val)), It(std::end(val))), [](char c) {
      return c == '\0';
    });
  }
  
  std::string encode64(const std::string &val) {
    using namespace boost::archive::iterators;
    using It = base64_from_binary<transform_width<std::string::const_iterator, 6, 8>>;
    auto tmp = std::string(It(std::begin(val)), It(std::end(val)));
    return tmp.append((3 - val.size() % 3) % 3, '=');
  }
  
  struct gnuplotImage
  {
    std::string imageText;
  };
  
  class GnuplotWidget
  {
  public:
    gnuplotImage plotData(const std::vector<double> &x,const std::vector<double> &y,std::istream* commands = nullptr);
    gnuplotImage plotData(const std::vector<double> &x,const std::vector<std::vector<double>> &y,std::istream* commands = nullptr);
    gnuplotImage plotData(const std::vector<std::vector<double>> &x,const std::vector<std::vector<double>> &y,std::istream* commands = nullptr);
    gnuplotImage plotData3D(const std::vector<double> &x,const std::vector<double> &y,const std::vector<std::vector<double>> &z,std::istream* commands = nullptr);
  };
  
}


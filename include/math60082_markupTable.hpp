#pragma once
#include <iomanip>
#include <iostream>

namespace MATH60082{
  
  template<class T=std::string>
  void addColumnEntry(const T& t=std::string("--------------"))
  {
    static const int columnWidth=14;
    std::cout << "|" << std::setw(columnWidth)  << t;
  }
  void endRow();
  void emptyTableRow(unsigned int n);
  
  template<typename T>
  void tableRow(T v) {
    addColumnEntry(v);
    endRow();
  }
  template<typename T, typename... Args>
  void tableRow(T first, Args... args) {
    addColumnEntry(first);
    tableRow(args...);
  }
  
}

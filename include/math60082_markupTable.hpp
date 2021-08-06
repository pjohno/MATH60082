#pragma once
#include <iomanip>
#include <iostream>
using namespace std;

namespace MATH60082{
    
  
  
  template<class T=std::string>
  void addColumnEntry(const T& t=std::string("--------------"))
  {
    static const int columnWidth=14;
    cout << "|" << setw(columnWidth)  << t;
  }
  void endRow()
  {
    cout << "|" <<endl;
  }
  void emptyTableRow(unsigned int n)
  {
    for(unsigned int i=0;i<n;i++)
    {
      addColumnEntry();
    }
    cout << "|" <<endl;
  }
  
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

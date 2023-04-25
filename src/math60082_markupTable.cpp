#include <iomanip>
#include <iostream>
#include "math60082_markupTable.hpp"
using namespace::std;

namespace MATH60082{
  
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
  
}

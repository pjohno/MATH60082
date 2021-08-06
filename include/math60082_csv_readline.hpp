#pragma once

#include <string>
#include <vector>
#include <iostream>

namespace MATH60082{
    
    /**
     * @brief Read a single row of comma seperated variables
     * @details This will return a vector of strings from single line read from a filestream
     * using stringstream to separate the columns.
     */
    struct ReadCSV
    {
      /// choice of delimiter
      static char delim;
      /** @brief Read a single row of comma seperated variables
       * @param input the input file stream
       * @return the vector of strings contained in each column of the row
     * @details An example usage would be
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
     * #include "mCSVReadline.hpp"
     * #include <iostream>
     * #include <fstream>
     * using namespace std;
     * 
     * int main()
     * {
     *   ifstream input("mCSVReadline-ex.csv");
     *   while(!(input.eof()))
     *   {
     *     vector<string> rowData = CommonLib::Math::ReadCSV::readRow(input);
     *     for(const auto& data : rowData)
     *     {
     *       cout << "|-" << data << "-";
     *     }
     *     cout << endl;
     *   }
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~
     * where the following input file
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~
     * the,cat, sat ,  on  ,the,, , ,   ,mat
     * , ,  ,   ,   
     * 1,1.2, 2.44 
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~
     * generates the output
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~
     * |-the-|-cat-|- sat -|-  on  -|-the-|--|- -|- -|-   -|-mat-
     * |--|- -|-  -|-   -|-   -
     * |-1-|-1.2-|- 2.44  -
     * 
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
      static std::vector<std::string> readRow(std::ifstream &input);
      
    };
    
  }


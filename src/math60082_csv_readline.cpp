#include "math60082_csv_readline.hpp"
#include <fstream>
#include <sstream>
using namespace MATH60082;

// #define DEBUG
char ReadCSV::delim=',';

std::vector<std::string> ReadCSV::readRow(std::ifstream &input)
{
    // create a new row of data
    std::vector<std::string> stringArray;
    // declare string buffer
    std::string buffer;
    //        get line from file
    std::getline(input, buffer, '\n');
    #ifdef DEBUG
    std::cout << " Buffer :: \"" << buffer << "\"\n";
    #endif
    if(buffer.size()==0)return stringArray;
    
    // read the buffer into a stringstream
    std::stringstream ss(buffer);
    while (ss)
    {
        std::string s(" ");
        if (!std::getline( ss, s, delim )) break;
        stringArray.push_back( s );
    }
    #ifdef DEBUG
    for(auto s : stringArray)
        std::cout << " s :: \"" << s << "\"\n";
    #endif
    return stringArray;  
}


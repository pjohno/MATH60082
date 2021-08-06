// include standard stuff
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
/**
 * @brief Read a single row of comma seperated variables
 * @details This will return a vector of strings from single line read from a filestream
 * using stringstream to separate the columns.
 */
struct ReadCSV
{
    /** @brief Read a single row of comma seperated variables
     * @param input the input file stream
     * @return the vector of strings contained in each column of the row
     */
    static std::vector<std::string> readRow(std::ifstream &input)
    {
        // create a new row of data
        std::vector<std::string> stringArray;
        // declare string buffer
        std::string buffer;
        //        get line from file
        std::getline(input, buffer, '\n');
        if(buffer.size()==0)return stringArray;
        
        // read the buffer into a stringstream
        std::stringstream ss(buffer);
        while (ss)
        {
            std::string s(" ");
            if (!std::getline( ss, s, ',' )) break;
            stringArray.push_back( s );
        }
        return stringArray;  
    }
};

int main()
{
    
    std::ofstream output("example.csv");
    for(int i=0;i<10;i++)
    {
        output << i << " , " << sin(i) << " , " << cos(i) << std::endl;
    }
    output.close();
    
    std::ifstream input("example.csv");
    while(!(input.eof()))
    {
        auto rowData = ReadCSV::readRow(input);
        for(const auto& data : rowData)
        {
            std::cout << "| " << data << " |";
        }
        std::cout << std::endl;
    }
}

#include "boost/process.hpp"
#include "math60082_gnuplot.hpp"
#include <iostream>
#include <iosfwd>
#include <boost/chrono.hpp>
#include <sstream>

namespace MATH60082
{
    void GnuplotWidget::plotToFile(std::istream* commands)
    {

        boost::process::opstream inputStream;
        boost::process::ipstream outputStream;
        boost::process::child myProcess("gnuplot",boost::process::std_out > outputStream,boost::process::std_in < inputStream);

        inputStream << "set terminal png"<< std::endl;
        if(commands){
            while(!commands->eof())
            {
                std::string line;
                getline(*commands,line);
                inputStream << line << std::endl;
            }
        }
        inputStream << "q"<< std::endl;
        myProcess.wait();
        std::string value;
        while(!outputStream.eof())
        {
            std::string line;
            getline(outputStream,line);
            value = value + line + '\n';
        }
        return;
    }

    gnuplotImage GnuplotWidget::plotCommand(std::istream* commands)
    {
        
        boost::process::opstream inputStream;
        boost::process::ipstream outputStream;
        boost::process::child myProcess("gnuplot",boost::process::std_out > outputStream,boost::process::std_in < inputStream);
        
        inputStream << "set terminal png"<< std::endl;
        if(commands){
            while(!commands->eof())
            {
                std::string line;
                getline(*commands,line);
                inputStream << line << std::endl;
            }
        } 
        inputStream << "q"<< std::endl;
        myProcess.wait();
        std::string value;
        while(!outputStream.eof())
        {
            std::string line;
            getline(outputStream,line);
            value = value + line + '\n';
        }
        return { value };
    }
    
    gnuplotGif GnuplotWidget::plotCommandGif(std::istream* commands)
    {


        boost::process::opstream inputStream;
        boost::process::ipstream outputStream;
        boost::process::child myProcess("gnuplot",boost::process::std_out > outputStream,boost::process::std_in < inputStream);

        inputStream << "set terminal gif"<< std::endl;
        if(commands){
            while(!commands->eof())
            {
                std::string line;
                getline(*commands,line);
                inputStream << line << std::endl;
            }
        }
        inputStream << "q"<< std::endl;
        myProcess.wait();
        std::string value;
        while(!outputStream.eof())
        {
            std::string line;
            getline(outputStream,line);
            value = value + line + '\n';
        }
        return { value };
    }

    gnuplotImage GnuplotWidget::plotData(const std::vector<double> &x,const std::vector<double> &y,std::istream* commands)
    {
        
        int xn = x.size();
        int yn = y.size();
        if(xn!=yn){
            std::cerr << " GnuplotWidget::plotData(x,y) := BAD SIZING OF VECTORS " << std::endl;
            return { "" };
        }
        
        boost::process::opstream inputStream;
        boost::process::ipstream outputStream;
        boost::process::child myProcess("gnuplot",boost::process::std_out > outputStream,boost::process::std_in < inputStream);
        
        inputStream << "set terminal png"<< std::endl;
        inputStream << "unset key"<< std::endl;
        if(commands){
            while(!commands->eof())
            {
                std::string line;
                getline(*commands,line);
                inputStream << line << std::endl;
            }
        }
        inputStream << "p '-' u 1:2 w l"<< std::endl;
        for(int i=0;i<xn;i++)
        {
            inputStream << x[i] << " " << y[i] << std::endl;
        }
        inputStream << "e"<< std::endl;
        inputStream << "q"<< std::endl;
        myProcess.wait();
        std::string value;
        while(!outputStream.eof())
        {
            std::string line;
            getline(outputStream,line);
            value = value + line + '\n';
        }
        return { value };
    }
    
    gnuplotImage GnuplotWidget::plotData(const std::vector<double> &x,const std::vector<std::vector<double>> &y,std::istream* commands)
    {
        
        int xn = x.size();
        for(int j=0;j<y.size();j++)
        {
            int yn = y[j].size();
            if(xn!=yn){
                std::cerr << " GnuplotWidget::plotData(x,y["<<j<<"]) := BAD SIZING OF VECTORS " << std::endl;
                return { "" };
            }
        }
        
        boost::process::opstream inputStream;
        boost::process::ipstream outputStream;
        boost::process::child myProcess("gnuplot",boost::process::std_out > outputStream,boost::process::std_in < inputStream);
        
        inputStream << "set terminal png"<< std::endl;
        inputStream << "unset key"<< std::endl;
        if(commands)
        {
            while(!commands->eof())
            {
                std::string line;
                getline(*commands,line);
                inputStream << line << std::endl;
            }
        }
        inputStream << "p";
        for(int j=0;j<y.size();j++)
        {
            if(j==0)
                inputStream << " '-' ";
            else
                inputStream << " '' ";
            inputStream << " u 1:2 w l";
            if(j<y.size()-1)
                inputStream << " , ";
        }
        inputStream  << std::endl;
        
        for(int j=0;j<y.size();j++)
        {
            for(int i=0;i<xn;i++)
            {
                inputStream << x[i] << " " << y[j][i] << std::endl;
            }
            inputStream << "e"<< std::endl;
        }
        inputStream << "q"<< std::endl;
        myProcess.wait();
        std::string value;
        while(!outputStream.eof())
        {
            std::string line;
            getline(outputStream,line);
            value = value + line + '\n';
        }
        return { value };
    }
    
    
    gnuplotImage GnuplotWidget::plotData(const std::vector<std::vector<double>> &x,const std::vector<std::vector<double>> &y,std::istream* commands)
    {
        int xTotalPlots = x.size();
        int yTotalPlots = y.size();
        if(xTotalPlots!=yTotalPlots){
                std::cerr << " GnuplotWidget::plotData(x,y) := BAD SIZING OF VECTORS " << std::endl;
                return { "" };
            }
        for(unsigned int j=0;j<y.size();j++)
        {
            int xn = x[j].size();
            int yn = y[j].size();
            if(xn!=yn){
                std::cerr << " GnuplotWidget::plotData(x["<<j<<"],y["<<j<<"]) := BAD SIZING OF VECTORS " << std::endl;
                return { "" };
            }
        }
        
        boost::process::opstream inputStream;
        boost::process::ipstream outputStream;
        boost::process::child myProcess("gnuplot",boost::process::std_out > outputStream,boost::process::std_in < inputStream);
        
        inputStream << "set terminal png"<< std::endl;
        inputStream << "unset key"<< std::endl;
        if(commands)
        {
            while(!commands->eof())
            {
                std::string line;
                getline(*commands,line);
                inputStream << line << std::endl;
            }
        }
        inputStream << "p";
        for(unsigned int j=0;j<y.size();j++)
        {
            if(j==0)
                inputStream << " '-' ";
            else
                inputStream << " '' ";
            inputStream << " u 1:2 w l";
            if(j<y.size()-1)
                inputStream << " , ";
        }
        inputStream  << std::endl;
        
        for(unsigned int j=0;j<y.size();j++)
        {
            for(unsigned int i=0;i<y[j].size();i++)
            {
                inputStream << x[j][i] << " " << y[j][i] << std::endl;
            }
            inputStream << "e"<< std::endl;
        }
        inputStream << "q"<< std::endl;
        myProcess.wait();
        std::string value;
        while(!outputStream.eof())
        {
            std::string line;
            getline(outputStream,line);
            value = value + line + '\n';
        }
        return { value };
    }
    
    
    gnuplotImage GnuplotWidget::plotData3D(const std::vector<double> &x,const std::vector<double> &y,const std::vector<std::vector<double>> &z,std::istream* commands)
    {
        
        int xn = x.size();
        int yn = y.size();
        int z_xn = z.size();
        if(xn!=z_xn){
            std::cerr << " GnuplotWidget::plotData(x,y) := BAD SIZING OF VECTORS " << std::endl;
            return { "" };
        }
        for(int i=0;i<z.size();i++)
        {
            int z_yn = z[i].size();
            if(yn!=z_yn){
                std::cerr << " GnuplotWidget::plotData(y,z["<<i<<"]) := BAD SIZING OF VECTORS " << std::endl;
                return { "" };
            }
        }
        
        boost::process::opstream inputStream;
        boost::process::ipstream outputStream;
        boost::process::child myProcess("gnuplot",boost::process::std_out > outputStream,boost::process::std_in < inputStream);
        
        inputStream << "set terminal png"<< std::endl;
        inputStream << "unset key"<< std::endl;
        if(commands){
            while(!commands->eof())
            {
                std::string line;
                getline(*commands,line);
                inputStream << line << std::endl;
            }
        }
        inputStream << "sp '-' u 1:2:3 w l"<< std::endl;
        for(int i=0;i<xn;i++)
        {
            for(int j=0;j<yn;j++)
            {
                inputStream << x[i] << " " << y[j] <<  " " << z[i][j] << std::endl;
            }
            inputStream << std::endl;
        }
        inputStream << "e"<< std::endl;
        inputStream << "q"<< std::endl;
        myProcess.wait();
        std::string value;
        while(!outputStream.eof())
        {
            std::string line;
            getline(outputStream,line);
            value = value + line + '\n';
        }
        return { value };
    }
    
}


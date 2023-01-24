#include "math60082_leastSquaresFit.hpp"
#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
using namespace std;

namespace MATH60082 
{
    
    void LeastSquaresFit::generateFit(int p,std::vector<DataPoint> &data)
    {
        int n=data.size();
        
        std::sort(data.begin(),data.end(),[](const DataPoint &a,const DataPoint &b){return (a.x<b.x);});
        xTotal=data.back().x+data.front().x;
        xRange=data.back().x-data.front().x;
        
        Eigen::MatrixXd A(n,p+1);
        Eigen::VectorXd b(n);
        for (int i = 0; i < n; i++)
        {
            double xi = data[i].x;
            xi = (2*xi-xTotal)/xRange;
            double yi = data[i].y;
            double Tn2 = 1.;
            A(i, 0)= Tn2;
            double Tn1 = xi;
            A(i, 1) = Tn1;
            for(int j=2;j<=p;j++)
            {
                double Tn = 2.*Tn1*xi - Tn2;
                A(i, j) = Tn;
                Tn2 = Tn1;
                Tn1 = Tn;
            }
            b(i)= yi;
        }
        
        Eigen::VectorXd x(p+1);
        
        x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
        
        P.init(p+1,0.5*(xTotal-xRange),0.5*(xTotal+xRange));
        for(int j=0;j<=p;j++)
        {
            P[j] = x(j);
        }
        
    }
    
    // print the fitting details
    ostream& operator<<(ostream& output,const LeastSquaresFit& fit)
    {
        output << " LeastSquaresFit of degree ";
        output << fit.P.size()-1 << " valid over the range x\\in["<<0.5*(fit.xTotal-fit.xRange) << ":"<< 0.5*(fit.xTotal+fit.xRange) <<"]";
        output << " :: { ";
        output << fit.P << " ";
        output << "} which gives a range of the function ";
        output << " f\\in["<<fit(0.5*(fit.xTotal-fit.xRange)) << ":"<< fit(0.5*(fit.xTotal+fit.xRange)) <<"]\n";
//         output << "# covariance matrix:\n";
//         for(unsigned int i=0;i<fit.C.size();i++)
//         {
//             for(unsigned int j=0;j<fit.C[i].size();j++)
//             {
//                 std::cout << " " << fit.C[i][j] << " " ;
//             }
//             std::cout << std::endl;
//         }
//         output << "# chisq = " <<  fit.chiSq << endl;
        return output;
    }
    
}

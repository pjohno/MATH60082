#include "math60082_fit.hpp"
#include "math60082_cheb.hpp"
#include "math60082_markupTable.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <Eigen/Dense>

namespace MATH60082
{
    int runFitEigen(std::vector<DataPoint> &data,int p)
    {
        int n=data.size();
        
        std::sort(data.begin(),data.end(),[](const DataPoint &a,const DataPoint &b){return (a.x<b.x);});
        double xTotal=data.back().x+data.front().x;
        double xRange=data.back().x-data.front().x;
        
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
        
        MATH60082::ChebyshevPolynomial P;
        P.init(p+1,0.5*(xTotal-xRange),0.5*(xTotal+xRange));
        for(int j=0;j<=p;j++)
        {
            P[j] = x(j);
        }
        
        std::cout << " Solution is " << P << std::endl;
        MATH60082::tableRow("x","y","P(x)");
        MATH60082::emptyTableRow(3);
        for(auto di : data)
        {
            MATH60082::tableRow(di.x,di.y,P(di.x));
        }

        return 0;
    }
}

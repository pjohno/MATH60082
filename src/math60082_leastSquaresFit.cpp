#include "math60082_leastSquaresFit.hpp"
#include <iostream>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_multifit.h>
using namespace std;

namespace MATH60082 
{
    
    void LeastSquaresFit::generateFit(int p,std::vector<DataPoint> &data)
    {
        // sort the data
        std::sort(data.begin(),data.end(),[](const DataPoint &a,const DataPoint &b){return (a.x<b.x);});
        xTotal=data.back().x+data.front().x;
        xRange=data.back().x-data.front().x;
        
        int i;
        gsl_matrix *X, *cov;
        gsl_vector *y, *w, *c;
        
        int n=data.size();
        
        X = gsl_matrix_alloc (n, p+1);
        y = gsl_vector_alloc (n);
        w = gsl_vector_alloc (n);
        
        c = gsl_vector_alloc (p+1);
        cov = gsl_matrix_alloc (p+1, p+1);
        
        std::sort(data.begin(),data.end(),[](const DataPoint &a,const DataPoint &b){return (a.x<b.x);});
        
        for (i = 0; i < n; i++)
        {
            double xi = data[i].x;
            xi = (2*xi-xTotal)/xRange;
            double yi = data[i].y;
            double wi = 1.;//data[i].w;
            double Tn2 = 1.;
            gsl_matrix_set (X, i, 0, Tn2);
            double Tn1 = xi;
            gsl_matrix_set (X, i, 1, Tn1);
            for(int j=2;j<=p;j++)
            {
                double Tn = 2.*Tn1*xi - Tn2;
                gsl_matrix_set (X, i, j, Tn);
                Tn2 = Tn1;
                Tn1 = Tn;
            }
            gsl_vector_set (y, i, yi);
            gsl_vector_set (w, i, wi);
        }
        
        {
            gsl_multifit_linear_workspace * work 
            = gsl_multifit_linear_alloc (n, p+1);
            gsl_multifit_wlinear (X, w, y, c, cov,
                                  &chiSq, work);
            gsl_multifit_linear_free (work);
        }
        
        P.init(p+1,0.5*(xTotal-xRange),0.5*(xTotal+xRange));
        for(int j=0;j<=p;j++)
        {
            P[j] = gsl_vector_get (c, j);
        }
        C.clear();
        C.resize(p+1,std::vector<double>(p+1));
        for(int i=0;i<=p;i++)
        {
            for(int j=0;j<=p;j++)
            {
                C[i][j] = gsl_matrix_get(cov, i, j);
            }
        }
        
        gsl_matrix_free (X);
        gsl_vector_free (y);
        gsl_vector_free (w);
        gsl_vector_free (c);
        gsl_matrix_free (cov);
        
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
        output << "# covariance matrix:\n";
        for(unsigned int i=0;i<fit.C.size();i++)
        {
            for(unsigned int j=0;j<fit.C[i].size();j++)
            {
                std::cout << " " << fit.C[i][j] << " " ;
            }
            std::cout << std::endl;
        }
        output << "# chisq = " <<  fit.chiSq << endl;
        return output;
    }
    
}

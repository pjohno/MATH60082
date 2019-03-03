#pragma once
#include "math60082_data.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_multifit.h>


namespace MATH60082
{
    
    int runFit(std::vector<DataPoint> &data,int p)
    {
        int i;
        double chisq;
        gsl_matrix *X, *cov;
        gsl_vector *y, *w, *c;
        
        int n=data.size();
        
        X = gsl_matrix_alloc (n, p);
        y = gsl_vector_alloc (n);
        w = gsl_vector_alloc (n);
        
        c = gsl_vector_alloc (p);
        cov = gsl_matrix_alloc (p, p);
        
        std::sort(data.begin(),data.end(),[](const DataPoint &a,const DataPoint &b){return (a.x<b.x);});
        
        for (i = 0; i < n; i++)
        {
            double xi = data[i].x;
            double yi = data[i].y;
            double wi = 1.;//data[i].w;
            double Tn2 = 1.;
            gsl_matrix_set (X, i, 0, Tn2);
            double Tn1 = xi;
            gsl_matrix_set (X, i, 1, Tn1);
            for(int j=2;j<p;j++)
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
            = gsl_multifit_linear_alloc (n, p);
            gsl_multifit_wlinear (X, w, y, c, cov,
                                  &chisq, work);
            gsl_multifit_linear_free (work);
        }
        
        #define C(i) (gsl_vector_get(c,(i)))
        #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
        
        {
            std::cout << " y = ";
            for(int j=0;j<p;j++)
            {
                std::cout << gsl_vector_get(c,(j))<<" T_"<<j<<"(x) ";
                if(j<p-1)std::cout << " + ";
            }
            std::cout << std::endl;
            printf ("# covariance matrix:\n");
            for(int i=0;i<p;i++)
            {
                for(int j=0;j<p;j++)
                {
                    std::cout << " " << gsl_matrix_get(cov,(i),(j)) << " " ;
                }
                std::cout << std::endl;
            }
            printf ("# chisq = %g\n", chisq);
        }
        
        std::ofstream output("test.dat");
        for (i = 0; i < n; i++)
        {
            double yi = data[i].y;
            output << " " << yi << " ";
            double xi = data[i].x;
            output << xi <<" ";
            double y=0.;
            {
                double Tn2 = 1.;
                y += gsl_vector_get(c,(0))*Tn2;
                double Tn1 = xi;
                y += gsl_vector_get(c,(1))*Tn1;
                double Tn = 2.*Tn1*xi - Tn2;
                y += gsl_vector_get(c,(2))*Tn;
            }
            output << y << std::endl;
        }
        gsl_matrix_free (X);
        gsl_vector_free (y);
        gsl_vector_free (w);
        gsl_vector_free (c);
        gsl_matrix_free (cov);
        
        return 0;
    }
}

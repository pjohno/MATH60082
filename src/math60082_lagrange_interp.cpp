#include "math60082_lagrange_interp.hpp"
#include <cmath>
using namespace std;


namespace MATH60082
{
    // A generic lagrange interpolation function
    double lagrangeInterpolation(const vector<double>& y,const vector<double>& x,double x0,unsigned int n)
    {
        if(x.size()<n)return lagrangeInterpolation(y,x,x0,x.size());
        if(n==0)throw;
        int nHalf = n/2;
        int jStar;
        double dx=x[1]-x[0];
        if(n%2==0)
            jStar = int((x0 - x[0])/dx) -(nHalf-1);
        else
            jStar = int((x0 - x[0])/dx+0.5)-(nHalf);
        jStar=std::max(0,jStar);
        jStar=std::min(int(x.size()-n),jStar);
        if(n==1)return y[jStar];
        double temp = 0.;
        for(unsigned int i=jStar;i<jStar+n;i++){
            double  int_temp;
            int_temp = y[i];
            for(unsigned int j=jStar;j<jStar+n;j++){
                if(j==i){continue;}
                int_temp *= ( x0 - x[j] )/( x[i] - x[j] );
            }
            temp += int_temp;
        }
        // end of interpolate
        return temp;
    }
}   

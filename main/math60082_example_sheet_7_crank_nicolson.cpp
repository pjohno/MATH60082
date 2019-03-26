#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;

int main()
{
    // declare and initialise Black Scholes parameters
    double S0=1.974,X=2.,T=1.,r=0.05,sigma=0.4;
    // declare and initialise grid paramaters 
    int iMax=4,jMax=4;
    // declare and initialise local variables (ds,dt)
    double S_max=2*X;
    double dS=S_max/jMax;
    double dt=T/iMax;
    // create storage for the stock price and option price (old and new)
    vector<double> S(jMax+1),vOld(jMax+1),vNew(jMax+1);
    // setup and initialise the stock price 
    for(int j=0;j<=jMax;j++)
    {
        S[j] = j*dS;
    }
    // setup and initialise the final conditions on the option price 
    for(int j=0;j<=jMax;j++)
    {
        vOld[j] = max(X-S[j],0.);
        vNew[j] = max(X-S[j],0.);
        cout << iMax << " " << j << " " << S[j] << " " << vNew[j] << " " << vOld[j] << endl;
    }
    // start looping through time levels
    for(int i=iMax-1;i>=0;i--)
    {
        // declare vectors for matrix equations
        vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
        // set up matrix equations a[j]=
        a[0] = 0.;b[0] = 1.;c[0] = 0.;
        d[0] = X*exp(-r*(iMax-i)*dt);
        for(int j=1;j<jMax;j++)
        {
            a[j] = 0.25*(sigma*sigma*j*j-r*j);
            b[j] =-0.5*sigma*sigma*j*j - 0.5*r - 1./dt;
            c[j] = 0.25*(sigma*sigma*j*j+r*j);
            d[j] =-a[j]*vOld[j-1] - (b[j]+2./dt)*vOld[j] - c[j]*vOld[j+1];
        }
        a[jMax] = 0.;b[jMax] = 1.;c[jMax] = 0.;
        d[jMax] = 0.;
        // solve with SOR method
        int sor,iterMax=10;
        double tol = 1.e-8,omega=1.;
        for(sor=0;sor<iterMax;sor++)
        {
            double error=0.;
            cout << sor << " == \n";
            // implement sor in here
            {
                double y = (d[0] - c[0]*vNew[1])/b[0];
                vNew[0] = vNew[0] + omega*(y-vNew[0]); 
                cout << " ( " << vNew[0] << " , ";
            }
            for(int j=1;j<jMax;j++)
            {
                double y = (d[j] - a[j]*vNew[j-1] - c[j]*vNew[j+1])/b[j];
                vNew[j] = vNew[j] + omega*(y-vNew[j]); 
                cout << vNew[j] << " , ";
            }
            {
                double y = (d[jMax] - a[jMax]*vNew[jMax-1])/b[jMax];
                vNew[jMax] = vNew[jMax] + omega*(y-vNew[jMax]); 
                cout << vNew[jMax] << " )\n";
            }
            // calculate residual norm ||r|| as sum of absolute values
            error += fabs(d[0] - b[0]*vNew[0] - c[0]*vNew[1]);
            for(int j=1;j<jMax;j++)
                error += fabs(d[j] - a[j]*vNew[j-1] - b[j]*vNew[j] - c[j]*vNew[j+1]);
            error += fabs(d[jMax] - a[jMax]*vNew[jMax-1] - b[jMax]*vNew[jMax]);
            // make an exit condition when solution found
            if(error<tol)
                break;
        }
        cout << " Solved after " << sor << " iterations.\n";
        // set old=new
        vOld=vNew;
    }
    // finish looping through time levels
    
    // output the estimated option price
    double optionValue;
    {
        int jStar=S0/dS;
        double sum=0.;
        sum+=(S0 - S[jStar])/dS * vNew[jStar+1];
        sum+=(S[jStar+1] - S0)/dS * vNew[jStar];
        optionValue = sum;
    }
    cout << "Value of the option V(S="<<S0<<") = " <<  optionValue << "\n";
}

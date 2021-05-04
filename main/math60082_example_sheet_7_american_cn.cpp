#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <chrono>
using namespace std;

double americanPut_CN(double S0,double X,double T,double r,double sigma,int iMax,int jMax,double SMax,double omega,double tol,int iterMax)
{
    double dS=SMax/jMax;
    double dt=T/iMax;
    // create storage for the stock price and option price (old and new)
    std::vector<double> S(jMax+1),vOld(jMax+1),vNew(jMax+1);
    // setup and initialise the stock price 
    for(int j=0;j<=jMax;j++)
    {
        S[j] = j*dS;
    }
    // reset initial condition
    for(int j=0;j<=jMax;j++)
    {
        vOld[j] = std::max(X-S[j],0.);
        vNew[j] = std::max(X-S[j],0.);
    }
    // run through timesteps
    for(int i=iMax-1;i>=0;i--)
    {
        // declare vectors for matrix equations
        std::vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
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
        int sor;
        for(sor=0;sor<iterMax;sor++)
        {
            double error=0.;
            // implement sor in here
            {
                double y = (d[0] - c[0]*vNew[1])/b[0];
                y = std::max(y,X-S[0]);
                y = vNew[0] + omega*(y-vNew[0]); 
                error+=(y-vNew[0])*(y-vNew[0]);
                vNew[0] = y;
            }
            for(int j=1;j<jMax;j++)
            {
                double y = (d[j] - a[j]*vNew[j-1] - c[j]*vNew[j+1])/b[j];
                y = std::max(y,X-S[j]);
                y = vNew[j] + omega*(y-vNew[j]); 
                error+=(y-vNew[j])*(y-vNew[j]);
                vNew[j] = y;
            }
            {
                double y = (d[jMax] - a[jMax]*vNew[jMax-1])/b[jMax];
                y = std::max(y,X-S[jMax]);
                y = vNew[jMax] + omega*(y-vNew[jMax]); 
                error+=(y-vNew[jMax])*(y-vNew[jMax]);
                vNew[jMax] = y;
            }    
            // make an exit condition when solution found
            if(error<tol*tol)
                break;
        }
        if(sor>=iterMax)
        {
            std::cout << " Error NOT converging within required iterations\n";
            std::cout.flush();
            throw;
        }
        // set old=new
        vOld=vNew;
    }
    int jStar=S0/dS;
    double sum=0.;
    sum+=(S0 - S[jStar])/dS * vNew[jStar+1];
    sum+=(S[jStar+1] - S0)/dS * vNew[jStar];
    return sum;
}

int main()
{
    // declare and initialise Black Scholes parameters
    double S0,X,T,r,sigma;
    // declare and initialise grid paramaters
    int iMax,jMax;
    // initialise Black Scholes parameters 
    S0=2.;X=2.;T=1.;r=0.05;sigma=0.4; 
    cout << " n    | V_cn     | error | CPU (time) | V_extrap \n"; 
    double valueOld=0.;
    double valueExtrapOld=0.;
    for(int k=1;k<=7;k++) 
    { 
        int n=10*pow(2,k); 
        iMax = n; 
        jMax = n; 
        
        auto start = chrono::system_clock::now(); 
        double value_cn = americanPut_CN(S0,X,T,r,sigma,iMax,jMax,5*X,1.5,1.e-10,10000); 
        auto end = chrono::system_clock::now(); 
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start); 
        double value_extrap = (4.*value_cn - valueOld)/3.; 
        double value_extrap2 = (4.*value_extrap - valueExtrapOld)/3.; 
        
        cout << n << " | " << value_cn<< " | "; 
        cout << value_cn - valueOld << " | " ; 
        cout << elapsed.count()/1000.<< " | " ; 
        cout << value_extrap << " | " << value_extrap2 << "\n"; 
        valueOld = value_cn;
        valueExtrapOld = value_extrap;
    }
}

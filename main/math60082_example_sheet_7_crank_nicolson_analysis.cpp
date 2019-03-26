#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

/* 
 * O *N INPUT:
 * S0          -- initial stock price
 * X           -- exercise (strike) price
 * T           -- Time to expiry (years)
 * r           -- interest rate (per annum)
 * sigma       -- volatility (per annum^12)
 * iMax        -- number of time steps
 * jMax        -- number of space steps
 * SMax        -- Maximum value of S
 * omega       -- is the relaxation parameter 
 * tol         -- is the tolerance level
 * iterMax     -- is maximum iterations
 * ON OUTPUT:
 * return      -- the value of a European put option at S=S0, t=0
 */
double europeanPut_CN(double S0,double X,double T,double r,double sigma,int iMax,int jMax,double SMax,double omega,double tol,int iterMax)
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
                vNew[0] = vNew[0] + omega*(y-vNew[0]); 
            }
            for(int j=1;j<jMax;j++)
            {
                double y = (d[j] - a[j]*vNew[j-1] - c[j]*vNew[j+1])/b[j];
                vNew[j] = vNew[j] + omega*(y-vNew[j]); 
            }
            {
                double y = (d[jMax] - a[jMax]*vNew[jMax-1])/b[jMax];
                vNew[jMax] = vNew[jMax] + omega*(y-vNew[jMax]); 
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

double normalDistribution(double x)
{
    return 0.5*erfc(-x/sqrt(2.));
}

// return the value of a put option using the black scholes formula
double europeanPut_exact(double S,double X,double T,double r,double sigma)
{
    if(fabs(T)<1.e-14) // check if we are at maturity
    {
        if(S<X)return X-S;
        else return 0;
    }
    if((T)<=-1.e-14)return 0.; // option expired
    if(X<1.e-14*S)return 0.; // check if strike << asset then exercise with certainty
    if(S<1.e-14*X)return X*exp(-r*(T)) - S; // check if asset << strike then worthless
    if(sigma*sigma*(T)<1.e-14) // check if variance very small then no diffusion
    {
        if(S<X*exp(-r*(T)))return X*exp(-r*(T)) - S;
        else return 0.;
    }
    // calculate option price
    double d1=(log(S/X) + (r+sigma*sigma/2.)*(T))/(sigma*sqrt(T));
    double d2=(log(S/X) + (r-sigma*sigma/2.)*(T))/(sigma*sqrt(T));
    return normalDistribution(-d2)*X*exp(-r*(T)) - normalDistribution(-d1)*S ;
}

int main()
{
    // declare and initialise Black Scholes parameters
    double S0,X,T,r,sigma;
    // declare and initialise grid paramaters 
    int iMax,jMax;
    double SMax;
    
    // initialise Black Scholes parameters
    S0=1.973;X=2.;T=1.;r=0.05;sigma=0.4;
    // initialise grid paramaters 
    iMax=40;jMax=40;SMax=2*X;
    
    cout << europeanPut_CN(S0,X,T,r,sigma,iMax,jMax,SMax,1,1.e-6,1000);
    cout << " " << europeanPut_exact(S0,X,T,r,sigma);
    
    cout << "\n\n Running with S_0="<<S0<<"\n";
    cout << " n    | V_cn     | V_exact  | error   \n";
    for(int k=1;k<=5;k++)
    {
        int n=10*pow(2,k);
        iMax = n;
        jMax = n;
        double value_exact = europeanPut_exact(S0,X,T,r,sigma);
        double value_cn = europeanPut_CN(S0,X,T,r,sigma,iMax,jMax,5*X,1.4,1.e-8,1000);
        cout << n << " | " << value_cn;
        cout << " | " << value_exact << " | " ;
        cout << value_cn - value_exact << "\n";
    }
    S0 = X;
    
    cout << "\n\n Running with S_0="<<S0<<"\n";
    cout << " n    | V_cn     | V_exact  | error   \n";
    for(int k=1;k<=5;k++)
    {
        int n=10*pow(2,k);
        iMax = n;
        jMax = n;
        double value_exact = europeanPut_exact(S0,X,T,r,sigma);
        double value_cn = europeanPut_CN(S0,X,T,r,sigma,iMax,jMax,5*X,1.4,1.e-8,1000);
        cout << n << " | " << value_cn;
        cout << " | " << value_exact << " | " ;
        cout << value_cn - value_exact << "\n";
    }
    
    cout << "\n\n Running with S_0="<<S0<<" and extrapolation.\n";
    double valueOld=0;
    double valueOlder=0;
    S0 = X;
    cout << " n    | V_cn     | V_exact  | error  | V_extrap | error \n";
    for(int k=1;k<=5;k++)
    {
        int n=10*pow(2,k);
        iMax = n;
        jMax = n;
        double value_exact = europeanPut_exact(S0,X,T,r,sigma);
        double value_cn = europeanPut_CN(S0,X,T,r,sigma,iMax,jMax,5*X,1.2,1.e-8,1000);
        double value_extrap = (4.*value_cn - valueOld)/3.;
        cout << n << " | " << value_cn;
        cout << " | " << value_exact << " | " ;
        cout << value_cn - value_exact << " | ";
        cout << value_extrap << " | ";
        cout << value_extrap - value_exact << "\n";
        valueOld = value_cn;
    }
}

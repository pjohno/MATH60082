#include "math60082_cn.hpp"

namespace MATH60082
{
    
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
}

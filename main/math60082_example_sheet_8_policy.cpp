#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <chrono>
using namespace std;

/* 
 *    ON INPUT:
 *    a, b and c -- are the diagonals of the matrix
 *    rhs        -- is the right hand side
 *    ON OUTPUT:
 *    a, b, c             -- unchanged
 *    rhs                 -- solution to Ax=b
 */
void thomasSolve(const std::vector<double> &a,const std::vector<double> &b_,const std::vector<double> &c,std::vector<double> &rhs)
{
    int n=a.size();
    std::vector<double> b(n);
    // initial first value of b
    b[0]=b_[0];
    for(int j=1;j<n;j++)
    {
        b[j]=b_[j]-c[j-1]*a[j]/b[j-1];  
        rhs[j]=rhs[j]-rhs[j-1]*a[j]/b[j-1];
    }
    // calculate solution
    rhs[n-1]=rhs[n-1]/b[n-1];
    for(int j=n-2;j>=0;j--)
        rhs[j]=(rhs[j]-c[j]*rhs[j+1])/b[j];
}

/* Template code for a Policy Iteration with Finite Difference
 *  American put option with
 *     V(S,T) = max( X-S , 0 )
 *  and policy dictates that either
 *      dV/dt + 1/2 sigma^2 S^2 d^2V/dS^2 + r S dV/dS - r V = 0
 *  or
 *         V = X - S 
 *  holds. Use SOR iteration to decide on policy at each space step.
 * 
 *  At S=0 we have
 *     V = X
 *  As S->infty, assume V->0.
 */
double AmericanPut_policy(double S0,double X,double T,double r,double sigma,int iMax,int jMax,double SMax,double tol,int iterMax)
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
        // solve with policy iterations
        int policy;
        for(policy=0;policy<iterMax;policy++)
        {
            // create new vectors containing the FD approximations
            std::vector<double> aHat(a),bHat(b),cHat(c),dHat(d);
            // decide on policy here
            {
                // get SOR value
                double y = (d[0] - c[0]*vNew[1])/b[0];
                // if SOR value suggests policy is to exercise, adjust matrix equations
                if( y < X - S[0] )
                {
                    aHat[0]=0.;bHat[0]=1.;cHat[0]=0.;dHat[0]=X - S[0];
                }
            }
            for(int j=1;j<jMax;j++)
            {
                // get SOR value
                double y = (d[j] - a[j]*vNew[j-1] - c[j]*vNew[j+1])/b[j];
                // if SOR value suggests policy is to exercise, adjust matrix equations
                if( y < X - S[j] )
                {
                    aHat[j]=0.;bHat[j]=1.;cHat[j]=0.;dHat[j]=X - S[j];
                }
            }
            {
                double y = (d[jMax] - a[jMax]*vNew[jMax-1])/b[jMax];
                // if SOR value suggests policy is to exercise, adjust matrix equations
                if( y < X - S[jMax] )
                {
                    aHat[jMax]=0.;bHat[jMax]=1.;cHat[jMax]=0.;dHat[jMax]=X - S[jMax];
                }
            }
            // now solve *Hat matrix problem with thomas
            thomasSolve(aHat,bHat,cHat,dHat);
            // dHat now contains next guess at solution
            // check for differences between dHat and vNew
            double error=0.;
            for(int j=1;j<jMax;j++)
                error += (dHat[j]-vNew[j])*(dHat[j]-vNew[j]);
            // update vNew
            vNew = dHat;
            // make an exit condition when policy is stable
            if(error<tol*tol)
                break;
        }
        if(policy>=iterMax)
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
    double S0 = 9.735, X = 10., T = 1., r = 0.05, sigma = 0.4;
    // declare and initialise grid paramaters 
    for(int n=10;n<=10000;n*=2)
    {
        int iMax = n, jMax = n;
        auto start = std::chrono::steady_clock::now(); 
        cout << n << " :: " << AmericanPut_policy(S0, X, T, r, sigma, iMax, jMax, 5.*X,1.e-8,10000);
        auto finish = std::chrono::steady_clock::now(); 
        auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
        cout << " :: ("<< elapsed.count()<< ")"<< endl;
    }
    /*
     * OUTPUT >> " n :: V :: Time (seconds) "
10 :: 1.21534 :: (5.0663e-05)
20 :: 1.42365 :: (3.4247e-05)
40 :: 1.47005 :: (0.000102486)
80 :: 1.47495 :: (0.000435619)
160 :: 1.4743 :: (0.00137424)
320 :: 1.47494 :: (0.00528259)
640 :: 1.475 :: (0.0213249)
1280 :: 1.475 :: (0.0856818)
2560 :: 1.47501 :: (0.417922)
5120 :: 1.47501 :: (2.00561)
     */
}

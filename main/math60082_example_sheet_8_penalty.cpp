#include <iostream>
#include <fstream>
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

/* Template code for the penalty method with Finite Difference
 *  American put option with
 *     V(S,T) = max( X-S , 0 )
 *  and
 *      dV/dt + 1/2 sigma^2 S^2 d^2V/dS^2 + r S dV/dS - r V  + rho max( (X-S) - P , 0) = 0
 *  with penalty rho->\infty. 
 * 
 *  At S=0 we have
 *     V = X
 *  As S->infty, assume V->0.
 */
double AmericanPut_penalty(double S0,double X,double T,double r,double sigma,int iMax,int jMax,double SMax,double rho,double tol,int iterMax)
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
        a[0] = 0.;b[0] = -1./dt - 0.5*r;c[0] = 0.;
        d[0] = (-1./dt + 0.5*r)*vOld[0];
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
        int penaltyIt;
        for(penaltyIt=0;penaltyIt<iterMax;penaltyIt++)
        {
            // create new vectors containing the FD approximations
            std::vector<double> aHat(a),bHat(b),cHat(c),dHat(d);
            // decide on policy here
            for(int j=0;j<=jMax;j++)
            {
                // if current value suggests apply penalty, adjust matrix equations
                if( vNew[j] < X - S[j] )
                {
                    bHat[j]=b[j]-rho;dHat[j]=d[j] - rho*(X - S[j]);
                }
            }
            // now solve *Hat matrix problem with thomas
            thomasSolve(aHat,bHat,cHat,dHat);
            // dHat now contains next guess at solution
            // check for differences between dHat and vNew
            double error=0.;
            for(int j=0;j<=jMax;j++)
                error += (dHat[j]-vNew[j])*(dHat[j]-vNew[j]);
            // update vNew
            vNew = dHat;
            // make an exit condition when solution is converged
            if(error<tol*tol)
            {
                // cout << " # solved after "<< penaltyIt << " iterations\n\n" << endl;
                break;
            }
        }
        if(penaltyIt>=iterMax)
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
        cout << n << " :: " << AmericanPut_penalty(S0, X, T, r, sigma, iMax, jMax, 5.*X,1.e8,1.e-8,10000);
        auto finish = std::chrono::steady_clock::now(); 
        auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
        cout << " :: ("<< elapsed.count()<< ")"<< endl;
    }
    /*
     * OUTPUT >> " n :: V :: Time (seconds) "
10 :: 1.21534 :: (5.5233e-05)
20 :: 1.42365 :: (4.5895e-05)
40 :: 1.47005 :: (0.000127086)
80 :: 1.47495 :: (0.000408222)
160 :: 1.4743 :: (0.00149963)
320 :: 1.47494 :: (0.00491252)
640 :: 1.475 :: (0.0199213)
1280 :: 1.475 :: (0.0794912)
2560 :: 1.47501 :: (0.390522)
5120 :: 1.47501 :: (1.89866)     
     */
}

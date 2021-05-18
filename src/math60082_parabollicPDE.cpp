
#include "math60082_parabollicPDE.hpp"
#include "math60082_tridag.hpp"
#include <algorithm>

namespace MATH60082
{
    
    int ParabollicPDE::solve(double theta,double xMin,double xMax,int jMax,double t,double T,int iMax,bool updateStencils)
    {
        // check grid setup correctly
        if(jMax<2){
            std::cout << " Grid not big enough to run finite difference " << std::endl;
            return 1;
        }
        // check grid setup correctly
        if(iMax<1){
            std::cout << " no time steps " << std::endl;
            return 2;
        }
        double dx=(xMax-xMin)/jMax;
        double dt=(T-t)/iMax;
        // create storage for the stock price and option price (old and new)
        x.resize(jMax+1);
        u.resize(jMax+1);
        // reset initial condition
        for(int j=0;j<=jMax;j++)
            x[j] = xMin + j*dx;
        for(int j=0;j<=jMax;j++)
            u[j] = IC(x[j]);
        #ifdef __MY_DEBUG__
        std::cout << " ## INITIAL CONDITION \n";
        for(int j=0;j<=jMax;j++)
            std::cout << x[j] << " " << u[j] << std::endl;
        #endif    
        
        VectorStencil L(jMax+1);
        
        double mTheta=1.-theta;
        // run through timesteps
        // declare vectors for matrix equations
        if(updateStencils){
            for(int i=0;i<iMax;i++)
            {
                L.createStencils(alpha,beta,x,(i+theta)*dt);
                std::vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
                FullStencil bcMin = BC_min(x[0],(i+1)*dt,u,dx,dt);
                a[0] = bcMin.a;b[0] = bcMin.b;c[0] = bcMin.c;d[0] = bcMin.d;
                for(int j=1;j<jMax;j++)
                {
                    a[j] = theta*L[j].a;
                    b[j] = -1./dt + theta*L[j].b + theta*gamma(x[j],(i+theta)*dt);
                    c[j] = theta*L[j].c;
                    d[j] = -mTheta*L[j].a*u[j-1] 
                    - ( 1./dt + mTheta*L[j].b + mTheta*gamma(x[j],(i+theta)*dt))*u[j] 
                    - mTheta*L[j].c*u[j+1];
                }
                FullStencil bcMax = BC_max(x[jMax],(i+1)*dt,u,dx,dt);
                a[jMax] = bcMax.a;b[jMax] = bcMax.b;c[jMax] = bcMax.c;d[jMax] = bcMax.d;
                // solve with SOR method
                thomasSolve(a,b,c,d);
                // set old=new
                u=d;
            }
        }
        else{
            L.createStencils(alpha,beta,x);
            std::vector<double> a(jMax+1),b(jMax+1),c(jMax+1);
            FullStencil bcMin = BC_min(x[0],0.,u,dx,dt);
                a[0] = bcMin.a;b[0] = bcMin.b;c[0] = bcMin.c;
                for(int j=1;j<jMax;j++)
                {
                    a[j] = theta*L[j].a;
                    b[j] = -1./dt + theta*L[j].b + theta*gamma(x[j],0);
                    c[j] = theta*L[j].c;
                }
                FullStencil bcMax = BC_max(x[jMax],0,u,dx,dt);
                a[jMax] = bcMax.a;b[jMax] = bcMax.b;c[jMax] = bcMax.c;
            for(int i=0;i<iMax;i++)
            {
                std::vector<double> d(jMax+1);
                FullStencil bcMin = BC_min(x[0],(i+1)*dt,u,dx,dt);
                d[0] = bcMin.d;
                for(int j=1;j<jMax;j++)
                {
                    d[j] = -mTheta*L[j].a*u[j-1] 
                    - ( 1./dt + mTheta*L[j].b + mTheta*gamma(x[j],(i+theta)*dt))*u[j] 
                    - mTheta*L[j].c*u[j+1];
                }
                FullStencil bcMax = BC_max(x[jMax],(i+1)*dt,u,dx,dt);
                d[jMax] = bcMax.d;
                // solve with SOR method
                thomasSolve(a,b,c,d);
                // set old=new
                u=d;
            }
            
        }
        #ifdef __MY_DEBUG__
        std::cout << " ## SOLUTION \n";
        for(int j=0;j<=jMax;j++)
            std::cout << x[j] << " " << u[j] << "\n";
        std::cout << std::endl;
        #endif   
        return 0;
    }
    
    std::ostream& operator<<(std::ostream& output,const ParabollicPDE& pde)
    {
        for(unsigned int j=0;j<pde.x.size();j++)
            output << pde.x[j] << " " << pde.u[j] << "\n";
        output << std::endl;
        return output;
    }
}

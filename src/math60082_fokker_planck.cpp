#include "math60082_fokker_planck.hpp"
#include "math60082_parabollicPDE.hpp"
#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>
#include <boost/math/distributions.hpp>
#include <chrono>
#include <fstream>
using namespace std;
namespace MATH60082{
    
    int runFokkerPlanck_Example1()
    {
        double x0=9.3574,xT=10.41265;
        double kappa=0.134,theta=10.1463,gam=0.3,T=1,sigma=1.5317;
        
        MATH60082::ParabollicPDE fockerPlank_OU;
        
        
        fockerPlank_OU.alpha = [&](double x,double t){return 2*sigma*sigma*gam*pow(x,2.*gam-1.)-kappa*(theta-x);};
        fockerPlank_OU.beta = [&](double x,double t){return sigma*pow(x,gam);};
        fockerPlank_OU.gamma = [&](double x,double t){return sigma*sigma*gam*(2*gam-1.)*pow(x,2.*gam-2.)+kappa;};
        
        fockerPlank_OU.IC = [&](double x){
            double dx=fockerPlank_OU.x[1]-fockerPlank_OU.x[0];
            int jStar = (x-fockerPlank_OU.x[0])/dx;
            int jStar0 = (x0-fockerPlank_OU.x[0])/dx;
            if(jStar==jStar0)
                return 1./dx;
            else
                return 0.;
        };
        fockerPlank_OU.BC_min =
        [&](double x,double t,const vector<double> &uOld,double dx,double dt){
            double xStar = x+0.5*dx;
            double flux_u=kappa*(theta-xStar) - sigma*sigma*gam*pow(xStar,2.*gam-1.);
            double flux_du=0.5*sigma*sigma*gam*pow(xStar,2.*gam);
            double a=0.;
            double b=flux_du/dx+0.5*flux_u;
            double c=-flux_du/dx+0.5*flux_u;
            double d=0.;
            return MATH60082::FullStencil{a,b,c,d};
        };
        fockerPlank_OU.BC_max =
        [&](double x,double t,const vector<double> &uOld,double dx,double dt){ double a=0.,b=1.,c=0.,d=0.;
            return MATH60082::FullStencil{a,b,c,d};
            
        };
        
        boost::math::normal OU(x0*exp(-kappa*T) + theta*(1-exp(-kappa*T)),sigma*sqrt((1-exp(-2.*kappa*T))/(2.*kappa)));
        std::cout.precision(8);
        std::cout << " Solution pdf(x=" << xT << ",T="<<T<<") = " << pdf(OU,xT) << endl;
        
        double valueOld=1.;
        for(int k=1;k<=10;k++)
        {
            int n=4*pow(2,k);
            auto start = std::chrono::steady_clock::now(); 
            double sd_est = sigma*pow(x0,gam)*sqrt((1-exp(-2.*kappa*T))/(2.*kappa));
            double xMin=max(0.,x0-10*sd_est);
            int error =  fockerPlank_OU.solve(0.5,xMin,x0+10*sd_est,2*n,0.,T,n,false);
            if(error>0)throw;
            
            auto finish = std::chrono::steady_clock::now(); 
            auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
            
            double value=fockerPlank_OU(xT);
            double valueExtrap=(4.*value - valueOld)/3.;
            std::cout << n << " " << value << " " ;
            std::cout << valueExtrap << " ";
            std::cout << " :: ("<< elapsed.count()<< ")"<<std::endl;
            valueOld=value;
            
            ofstream output("test.dat");
            output << fockerPlank_OU;
            
        }
        
        return 0;
    }
}

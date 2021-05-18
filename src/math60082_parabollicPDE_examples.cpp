#include "math60082_parabollicPDE.hpp"
#include "math60082_parabollicPDE_examples.hpp"
#include "math60082_blackScholes.hpp"
#include "math60082_vasicekModel.hpp"
#include "doxygen_table.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <chrono>
#include <boost/math/distributions.hpp>
#include <boost/math/tools/toms748_solve.hpp>
using namespace std;

namespace MATH60082{
    
    int runParabollicPDE_Example1()
    {
        double S0=10;
        double X=10.1;
        double r=0.05,D=0.02,T=1,sigma=0.2;
        
        MATH60082::ParabollicPDE fockerPlank_geometric;
        
        fockerPlank_geometric.alpha = [&](double S,double t){return (2*sigma*sigma-r+D)*S;};
        fockerPlank_geometric.beta = [&](double S,double t){return sigma*S;};
        fockerPlank_geometric.gamma = [&](double S,double t){return sigma*sigma-r+D;};
        fockerPlank_geometric.IC = [&](double S){
            double dx=fockerPlank_geometric.x[1]-fockerPlank_geometric.x[0];
            int jStar = (S-fockerPlank_geometric.x[0])/dx;
            int jStar0 = (S0-fockerPlank_geometric.x[0])/dx;
            if(jStar==jStar0)
                return 1./dx;
            else
                return 0.;
        };
        fockerPlank_geometric.BC_min =
        [&](double S,double t,const vector<double> &uOld,double dS,double dt){
            double a=0.,b=1.,c=0.,d=0.;
            return MATH60082::FullStencil{a,b,c,d};};
            fockerPlank_geometric.BC_max =
            [&](double S,double t,const vector<double> &uOld,double dS,double dt){
                double a=0.,b=1.,c=0.,d=0.;
                return MATH60082::FullStencil{a,b,c,d};};
                
                boost::math::lognormal LN(log(S0) + (r-D-0.5*sigma*sigma)*T,sigma*sqrt(T));
                std::cout.precision(8);
                std::cout << " Solution pdf(S=" << X << ",T="<<T<<") = " << pdf(LN,X) << endl;
                
                double valueOld=1.;
                for(int k=1;k<=8;k++)
                {
                    int n=4*pow(2,k);
                    auto start = std::chrono::steady_clock::now(); 
                    int error =  fockerPlank_geometric.solve(0.5,0.,5*S0,5*n,0.,T,n,false);
                    if(error>0)throw;
                    
                    auto finish = std::chrono::steady_clock::now(); 
                    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
                    
                    double value=fockerPlank_geometric(X);
                    double valueExtrap=(4.*value - valueOld)/3.;
                    std::cout << n << " " << value << " " ;
                    std::cout << valueExtrap << " ";
                    std::cout << " :: ("<< elapsed.count()<< ")"<<std::endl;
                    valueOld=value;
                }
                
                return 0;
    }
    
    int runParabollicPDE_Example2()
    {
        double S0=9.735;
        double r=0.05,D=0.02,X=10.,T=1.,sigma=0.3;
        MATH60082::ParabollicPDE blackScholesCall;
        blackScholesCall.alpha= [&](double S,double t){return (r-D)*S;};
        blackScholesCall.beta=     [&](double S,double t){return sigma*S;};
        blackScholesCall.gamma=     [&](double S,double t){return -r;};
        blackScholesCall.IC=     [&](double S){return max(S-X,0.);};
        blackScholesCall.BC_min=     [&](double S,double t,const vector<double> &uOld,double dS,double dt){
            double a=0.,b=1.,c=0.,d=0.;
            return MATH60082::FullStencil{a,b,c,d};};
            blackScholesCall.BC_max=         [&](double S,double t,const vector<double> &uOld,double dS,double dt){
                double a=0.,b=1.,c=0.,d=S*exp(-D*(T-t)) - X*exp(-r*(T-t));
                return MATH60082::FullStencil{a,b,c,d};};
                
                std::cout.precision(8);
                std::cout << " Solution C(S=" << S0 << ",t="<<0<<";X="<<X<<",r="<<r<<",D="<<D<<",sigma="<<sigma<<") = " << MATH60082::blackScholesCallOption(S0,X,T,r,D,sigma) << endl;
                
                double valueOld=1.;
                for(int k=1;k<=8;k++)
                {
                    int n=4*pow(2,k);
                    auto start = std::chrono::steady_clock::now(); 
                    int error =  blackScholesCall.solve(0.5,0.,5*X,5*n,0.,T,n,false);
                    if(error>0)throw;
                    auto finish = std::chrono::steady_clock::now(); 
                    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
                    double value=blackScholesCall(S0);
                    std::cout << n << " | " << value << " | " ;
                    std::cout << (4.*value - valueOld)/3. << " | ";
                    std::cout << " | ("<< elapsed.count()<< ") |"<<std::endl;
                    valueOld=value;
                }
                
                return 0;
    }
    
    int runParabollicPDE_Example3()
    {
        double r0=0.0135;
        double kappa=0.05,theta=0.02,T=1.,sigma=0.07;
        MATH60082::ParabollicPDE VasicekBond;
        
        VasicekBond.alpha = [&](double r,double t){return kappa*(theta-r);};
        VasicekBond.beta =     [&](double r,double t){return sigma;};
        VasicekBond.gamma =     [&](double r,double t){return -r;};
        VasicekBond.IC =     [&](double r){return 1.;};
        VasicekBond.BC_min =    [&](double r,double t,const vector<double> &uOld,double dr,double dt){
            double a=0.;
            double b=1./dt+kappa*(theta-r)/dr+0.5*r;
            double c=-kappa*(theta-r)/dr;
            double d=(1./dt-kappa*(theta-r)/dr-0.5*r)*uOld[0] + (kappa*(theta-r)/dr)*uOld[1];
            return MATH60082::FullStencil{a,b,c,d};
        },
        VasicekBond.BC_max =     [&](double r,double t,const vector<double> &uOld,double dr,double dt){
            int jMax=uOld.size()-1;
            double a=-kappa*(theta-r)/dr;
            double b=1./dt+kappa*(theta-r)/dr+0.5*r;
            double c=0.;
            double d=(1./dt-kappa*(theta-r)/dr-0.5*r)*uOld[jMax] + (kappa*(theta-r)/dr)*uOld[jMax-1];
            return MATH60082::FullStencil{a,b,c,d};
        };
        
        std::cout.precision(8);
        MATH60082::VasicekInterestRateModel VM = VasicekInterestRateParameters{ kappa,theta,sigma};
        std::cout << " Solution B(r=" << r0 << ",t="<<0<<";kappa="<<kappa<<",theta="<<theta<<",sigma="<<sigma<<") = " << VM.discountBond(r0,0,T) << endl;
        
        double valueOld=1.;
        for(int k=1;k<=8;k++)
        {
            int n=4*pow(2,k);
            auto start = std::chrono::steady_clock::now(); 
            int error =  VasicekBond.solve(0.5,r0-5*T*sigma,r0+5*T*sigma,5*n,0.,T,n,false);
            if(error>0)throw;
            auto finish = std::chrono::steady_clock::now(); 
            auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
            double value=VasicekBond(r0);
            std::cout << n << " | " << value << " | " ;
            std::cout << (4.*value - valueOld)/3. << " | ";
            std::cout << " | ("<< elapsed.count()<< ") |"<<std::endl;
            valueOld=value;
        }
        
        double T2=1.,X=1.;
        MATH60082::ParabollicPDE BondOption;
        BondOption.alpha= [&](double r,double t){return kappa*(theta-r);};
        BondOption.beta=     [&](double r,double t){return sigma;};
        BondOption.gamma=     [&](double r,double t){return -r;};
        BondOption.IC=     [&](double r){return std::max(VasicekBond(r)-X,0.);};
        BondOption.BC_min=     [&](double r,double t,const vector<double> &uOld,double dr,double dt){
            double a=0.;
            double b=1./dt+kappa*(theta-r)/dr+0.5*r;
            double c=-kappa*(theta-r)/dr;
            double d=(1./dt-kappa*(theta-r)/dr-0.5*r)*uOld[0] + (kappa*(theta-r)/dr)*uOld[1];
            return MATH60082::FullStencil{a,b,c,d};
        };
        BondOption.BC_max=     [&](double r,double t,const vector<double> &uOld,double dr,double dt){
            int jMax=uOld.size()-1;
            double a=-kappa*(theta-r)/dr;
            double b=1./dt+kappa*(theta-r)/dr+0.5*r;
            double c=0.;
            double d=(1./dt-kappa*(theta-r)/dr-0.5*r)*uOld[jMax] + (kappa*(theta-r)/dr)*uOld[jMax-1];
            return MATH60082::FullStencil{a,b,c,d};
        };
        
        double rStar=r0;
        {
            double tol=1.-8;
            boost::uintmax_t iters=100;
            std::pair<double,double> sol = boost::math::tools::toms748_solve(
                [&](double r){return VasicekBond(r)-X;},r0-5*T*sigma,r0+5*T*sigma,[tol](const double& a, const double& b){return (fabs(a)<tol) && (fabs(a)<tol);},iters);
            rStar=0.5*(sol.first+sol.second);
        }
        
        for(int k=1;k<=8;k++)
        {
            int n=4*pow(2,k);
            auto start = std::chrono::steady_clock::now(); 
            int error =  BondOption.solve(0.5,rStar-5*T*sigma,rStar+5*T*sigma,5*n,0.,T,n,false);
            if(error>0)throw;
            auto finish = std::chrono::steady_clock::now(); 
            auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
            double value=BondOption(r0);
            std::cout << n << " | " << value << " | " ;
            std::cout << (4.*value - valueOld)/3. << " | ";
            std::cout << " | ("<< elapsed.count()<< ") |"<<std::endl;
            valueOld=value;
        }
        
        return 0;
    }
    
}

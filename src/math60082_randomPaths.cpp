#include "math60082_randomPaths.hpp"

namespace MATH60082
{
    
    std::vector< std::vector<double> > generateRandomPaths(
        int timesteps,double S0,double dt,double mu,double sigma,int N
    )
    {
        static std::mt19937 rng;
        double A = (mu - 0.5*sigma*sigma)*dt;
        double B = sigma*sqrt(dt);
        std::normal_distribution<double> Phi(A,B);
        std::vector<std::vector<double>> paths(N,std::vector<double> (timesteps+1));
        for(int i=0;i<N;i++)
        {
            paths[i][0] = S0;
            for(int k=1;k<=timesteps;k++)
            {
                paths[i][k] = paths[i][k-1] * exp(Phi(rng)); 
            }
        }
        return paths;
    }
    
    std::vector< std::vector<double> > generateLongstaffPaths(
        int &timesteps,double &dt,int &N
    )
    {
        timesteps=3;
        dt=1.;
        N=8;
        std::vector<std::vector<double>>paths = {
            { 1. , 1.09 , 1.08 , 1.34 } ,
            { 1. , 1.16 , 1.26 , 1.54 } ,
            { 1. , 1.22 , 1.07 , 1.03 } ,
            { 1. , 0.93 , 0.97 , 0.92 } ,
            { 1. , 1.11 , 1.56 , 1.52 } ,
            { 1. , 0.76 , 0.77 , 0.90 } ,
            { 1. , 0.92 , 0.84 , 1.01 } ,
            { 1. , 0.88 , 1.22 , 1.34 } ,
        };
        return paths;
    }
    
}

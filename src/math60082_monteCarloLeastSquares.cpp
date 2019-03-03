#include "math60082_randomPaths.hpp"
#include "math60082_data.hpp"
#include "math60082_leastSquaresFit.hpp"
#include "math60082_monteCarloLeastSquares.hpp"

namespace MATH60082
{
    
    double monteCarloLeastSquares(double S0,double maturity,double strikePrice,
                                  double interestRate,double sigma,int timeSteps,int M_degree,int n)
    {
        
        // local parameters
        double dt=maturity/timeSteps;
        
        std::vector<std::vector<double>> paths=generateRandomPaths(timeSteps,S0,dt,interestRate,sigma,n);
        std::vector<int> exercisePeriod(n,timeSteps);
        
        // assume that time period is t=2
        for(int k=timeSteps-1;k>=1;k--)
        { 
            // create an empty set of data points
            std::vector<DataPoint> data;
            data.reserve(paths.size());
            // go through each path, check if it is in the money
            // if it is add it as a datapoint with appropriate payoff
            for(int i=0;i<int(paths.size());i++)
            {
                if(paths[i][k]<strikePrice)
                {
                    int kStar=exercisePeriod[i];
                    double tStar=(kStar-k)*dt;
                    data.push_back({paths[i][k],std::max(strikePrice - paths[i][kStar],0.)*exp(-interestRate*tStar),1.});
                }
            }
            
            // get the continuation function
            LeastSquaresFit V;      
            V.generateFit(M_degree,data);
            
            // then select those paths and check against continuation function
            for(int i=0;i<int(paths.size());i++)
            {
                if(paths[i][k]<strikePrice)
                {
                    double exerciseValue=std::max(strikePrice - paths[i][k],0.);
                    if(exerciseValue > V(paths[i][k]))
                    {
                        exercisePeriod[i]=k;
                    }
                }
            }
            
        }
        
        // implement Monte Carlo on the paths
        double sum=0.;
        for(int i=0;i<n;i++)
        {
            // get early exercise time for American option
            int k=exercisePeriod[i];
            double exerciseValue=std::max(strikePrice - paths[i][k],0.);
            double exerciseTime=k*dt;
            sum+=exerciseValue*exp(-interestRate*exerciseTime);
        }
        return sum/n;
        
    }
    
}

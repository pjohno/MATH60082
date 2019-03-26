#include "math60082_leastSquaresFit.hpp"
#include <random>
using namespace MATH60082;
using namespace std;

int main()
{
    
    // number of exercise dates
    int timeSteps=3;
    // some parameters in the problem
    double maturity=3.,interestRate=0.06,strikePrice=1.1;
    // local parameters
    double dt=maturity/timeSteps;
    
    // the paths as given in Longstaff & Schartz (2001)
    vector<vector<double>> paths = 
    {
        { 1. , 1.09 , 1.08 , 1.34 } ,
        { 1. , 1.16 , 1.26 , 1.54 } ,
        { 1. , 1.22 , 1.07 , 1.03 } ,
        { 1. , 0.93 , 0.97 , 0.92 } ,
        { 1. , 1.11 , 1.56 , 1.52 } ,
        { 1. , 0.76 , 0.77 , 0.90 } ,
        { 1. , 0.92 , 0.84 , 1.01 } ,
        { 1. , 0.88 , 1.22 , 1.34 } ,
    };
    
    // all paths initially exected to be exercised at expiry
    vector<int> exercisePeriod = { 3 , 3 , 3 , 3 , 3 , 3 , 3 , 3 };
    
    // assume that time period is t=2
    for(int k=2;k>=1;k--)
    { 
        // create an empty set of data points
        vector<DataPoint> data;
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
        
        cout << "# Points for fitting\n x_i  y_i \n";
        for(auto di : data)
            cout << di.x << " " << di.y << endl;
        
        // get the continuation function
        LeastSquaresFit V;
        
        V.generateFit(2,data);
        
        cout << V << endl;
        
        cout << "# fitted estimation \n  x_i \t   y_i  \t V(x_i) \n";
        
        for(auto di : data)
            cout << di.x << "\t " << di.y << "\t " << V(di.x) << endl;
        
        cout << "# Cash flow at "<<k<<"th step: \n kstar  V \n";
        // then select those paths and check against continuation function
        for(int i=0;i<int(paths.size());i++)
        {
            
            if(paths[i][k]<strikePrice)
            {
                double exerciseValue=std::max(strikePrice - paths[i][k],0.);
                if(exerciseValue>V(paths[i][k]))
                {
                    exercisePeriod[i]=k;
                    cout << exercisePeriod[i] << " " << strikePrice - paths[i][k] << endl;
                }
                else 
                {
                    int kStar=exercisePeriod[i];
                    double tStar=(kStar-k)*dt;
                    cout << exercisePeriod[i] << " " << std::max(strikePrice - paths[i][kStar],0.)*exp(-interestRate*tStar) << endl;
                } 
            }
            else
            {
                int kStar=exercisePeriod[i];
                double tStar=(kStar-k)*dt;
                cout << exercisePeriod[i] << " " << std::max(strikePrice - paths[i][kStar],0.)*exp(-interestRate*tStar) << endl;
            }
        }
        
    }
    
    int n=paths.size();
    // implement Monte Carlo on the paths
    double sum=0.;
    for(int i=0;i<n;i++)
    {
        int k=exercisePeriod[i];
        double exerciseValue=std::max(strikePrice - paths[i][k],0.);
        double exerciseTime=k*dt;
        sum+=exerciseValue*exp(-interestRate*exerciseTime);
    }
    cout << " Option Value is " << sum/n << endl;
}



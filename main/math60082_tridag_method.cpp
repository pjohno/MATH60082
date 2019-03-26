#include <iostream>
#include <cmath>
#include <vector>
#include "math60082_sor.hpp"
#include "math60082_tridag.hpp"
using namespace std;

int main()
{
    // declare vector to store matrix A
    vector<double> a = { 0.,-1,-1,-1};
    vector<double> b = { 1.,2,2,1};
    vector<double> c = { 0.,-1,-1,0.};
    // declare vector for the rhs v and solution x 
    vector<double> rhs = {1.,0.25,0.5,0.};
    vector<double> x;
    // first try
    // reset initial guess for x
    x = {0.,0.,0.,0.};
    int sor;
    // iterMax=50, tol=1.e-6, omega = 1
    MATH60082::sorSolve(a,b,c,rhs,x,50,1.e-6,1.,sor);
    cout << "Converged after " << sor+1 << " iterations.\n";
    for(int varyOmega=5;varyOmega<=20;varyOmega++)
    {
        // reset initial guess for x
        x = {0.,0.,0.,0.};
        // this time omega = 1.1
        double omega = varyOmega/10.;
        MATH60082::sorSolve(a,b,c,rhs,x,5000,1.e-6,omega,sor);
        cout << "Try with omega = " << omega << " we get x = { ";
        for(auto xi : x)cout << xi << " ";
        cout << "} after " << sor+1 << " iterations.\n";
    }
    
    // reset x with the rhs
    x = rhs;
    // on exit x will contain the solution
    MATH60082::thomasSolve(a,b,c,x);
    cout << "Try with thomas solver we get x = { ";
        for(auto xi : x)cout << xi << " ";
        cout << "} " << endl;

}


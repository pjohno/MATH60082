#include <iostream>
#include <cmath>
#include <vector>
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

int main()
{
    // declare vector to store matrix A
    vector<double> a = { 0.,-1,-1,-1};
    vector<double> b = { 1.,2,2,1};
    vector<double> c = { 0.,-1,-1,0.};
    // declare vector for the rhs v and solution x 
    vector<double> rhs = {1.,0.25,0.5,0.};
    vector<double> x;
    // first reset x with the rhs
    x = rhs;
    // on exit x will contain the solution
    thomasSolve(a,b,c,x);
    cout << "Solve with thomas algorithm we get x = { ";
    for(auto xi : x)cout << xi << " ";
    cout << "} " << endl;
    
}


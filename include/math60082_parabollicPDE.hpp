
#pragma once
#include <iostream>
#include <vector>
#include <functional>
#include "math60082_stencils.hpp"
#include "math60082_lagrange_interp.hpp"

namespace MATH60082{
    
    /** @brief Structure to pass functions to parabollic solver algorithm
     * @details Solve a problem of the form
     * \f[ \frac{\partial u}{\partial t} = \alpha(x,t)\frac{\partial u}{\partial x} + \frac12\beta^2(x,t) \frac{\partial^2 u}{\partial x^2} + \gamma(x,t)u   \f]
     * where 
        **/
    struct ParabollicPDE
    {
        /** vector for the grid variable \f$ x \f$  **/
        std::vector<double> x;
        /** vector for the solution \f$ u \f$  **/
        std::vector<double> u;
        /** Here \f$\alpha\f$ is a function of \f$x\f$ and \f$t\f$ multiplying the first derivative in \f$x\f$. **/
        std::function<double(double,double)> alpha;
        /** Here \f$\frac12\beta^2\f$ is a function of \f$x\f$ and \f$t\f$ multiplying the second derivative in \f$x\f$. **/
        std::function<double(double,double)> beta;
        /** Here \f$\gamma\f$ is a function of \f$x\f$ and \f$t\f$ multiplying \f$u\f$. **/
        std::function<double(double,double)> gamma;
        /** The initial condition in the problem
         * \f[
         * u(x,t=0) = f(x) .
         * \f]
         * Given as function of \f$x\f$.
         **/
        std::function<double(double)> IC;
        /** The boundary condition at \f$x=x_\text{min}\f$, given as a stencil \f$\{a,b,c,d\}\f$ such that
         * \f[
         * b u^{i+1}_0 + c u^{i+1}_1 = d.
         * \f]
         * Given as function of \f$x\f$.
         * 
         **/
        std::function<FullStencil(double,double,const std::vector<double> &,double,double)> BC_min;
        /**
         * 
         **/
        std::function<FullStencil(double,double,const std::vector<double> &,double,double)> BC_max;
        /** @brief Return value of the function \f$ u(x=x_0,t=T) \f$.
        * @details Carries out a Lagrange interpolation at \f$ u(x) \f$, given the value of the function is stored as \f$ \mathbf{u}=\{ u_0,u_1,\dots, u_\text{jMax} \} \f$ at the grid points \f$ \mathbf{x}=\{ x_0,x_1,\dots, x_\text{jMax} \} \f$. The algorithm assumes an equally spaced grid is used to locate the nearest points in the grid.
        * @param x0 the value in the \f$ x \f$ grid to return a value
        * @param n number of points over which to interpolate
         **/
        double operator()(double x0,unsigned int n=4) const {return lagrangeInterpolation(u,x,x0,n);}
        /**
         * 
         * @param theta the scheme parameter, with \f$ \theta=0 \f$ corresponding to explicit, \f$ \theta=0 \f$ implicit and  \f$ \theta=\frac12 \f$ a Crank-Nicolson scheme.
         * @param xMin the location of the minimum x grid point
         * @param xMax the location of the maximum x grid point, should satisfy \f$ xMax > xMin \f$
         * @param jMax the number of divisions in the grid. There are  \f$ jMax+1 \f$ grid points, and is an integer \f$ jMax>2 \f$.
         * @param t the initial time value.
         * @param T the final time value.
         * @param iMax the number of time steps. Is an integer \f$ iMax\geq 1 \f$.
         * @param updateStencils update stencils at each time step. If stencils are constant in time, set this to false to speed up calculation.
         **/
        int solve(double theta,double xMin,double xMax,int jMax,double t,double T,int iMax,bool updateStencils=true);
        /** @brief Output the solution
        * @details 
        * @param output stream to send results to
        * @param pde the object we wish to print
        **/
        friend std::ostream& operator<<(std::ostream& output,const ParabollicPDE& pde);
        
        ParabollicPDE():x(100000),u(100000),alpha(nullptr),beta(nullptr),gamma(nullptr),IC(nullptr),BC_min(nullptr),BC_max(nullptr){}
        ~ParabollicPDE(){}
    };
    
}

#pragma once


namespace MATH60082{
    /** @brief Examples solving linear PDEs.
     * @details Solve a problem of the formFor Fokker Planck equation on SDE
     */
    /** @addtogroup Examples
     * @{
     */
    
    
    /**
     * Assume that we have an SDE of the form
     * \f[
     * dX = \kappa ( \theta - X ) + \sigma X^\gamma dW
     * \f]
     * and we wish to find the pdf function \f$u(x_T,T| x_t,t)\f$ which is the probability density of being at \f$x_T\f$ at time \f$T\f$ given that you started at \f$x_t\f$ at time \f$x_t\f$. It can be shown that \f$f\f$ satisfies the PDE
     * \f[
     * \frac{\partial u}{\partial t} + \frac{\partial}{\partial x}\{ \kappa(\theta - x) u \} 
     * -  \frac{\partial^2}{\partial x^2} \left\{ \frac12 \sigma^2 x^{2\gamma} u \right\} = 0
     * \f]
     * with initial condition
     * \f[
     * u(x,t=0) = \delta(x-x_0)
     * \f]
     * at \f$t=0\f$.
     * 
     * Note that the PDE can be rewritten as
     * \f[
     * \frac{\partial u}{\partial t} + \frac{\partial}{\partial x}\left\{ \kappa(\theta - x) u - \frac12 \sigma^2\frac{\partial }{\partial x} \{x^{2\gamma} u \}\right\} = 0,
     * \f]
     * and further simplified to 
     * \f[
     * \frac{\partial u}{\partial t} = \left( \kappa(\theta - x) - \sigma^22\gamma x^{2\gamma-1}\right)\frac{\partial u}{\partial x} + \frac12 \sigma^2 x^{2\gamma} \frac{\partial^2 u }{\partial x^2} + \left( \kappa + \frac12\sigma^2 2\gamma(2\gamma-1)x^{2\gamma-2}  \right) u.
     * \f]
     * We choose to solve the latter so that the coefficients from the resulting finite difference stencils are much simpler. This makes it easier to construct a stable scheme.
     * 
     * 
     * For the boundary condition on small \f$x\f$, we choose to put a flux condition in which means we solve
     * \f[
     * \frac12 \sigma^2\frac{\partial u}{\partial x} + \left( \sigma^2\gamma x^{2\gamma-1} - \kappa(\theta - x) \right)u
     * \f]
     */
    int runFokkerPlanck_Example1();
    
    
    /**  @}   */
    
} 

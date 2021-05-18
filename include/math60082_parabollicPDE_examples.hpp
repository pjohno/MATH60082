#pragma once


namespace MATH60082{
    
    /** @brief Examples solving linear PDEs.
     * @details Solve a problem of the form
     * \f[
\frac{\partial u}{\partial t} =
\alpha(x,t)\frac{\partial u}{\partial x}
+ \frac12\beta^2(x,t)\frac{\partial^2 u}{\partial x^2}
\f]
     */
    /** @addtogroup Examples
     * @{
     */
    
    /** @brief Example solving Fokker Planck for Geomtric Brownian Motion, using ParabollicPDE.
     * @details Solve Fokker Planck equation for the Geometric SDE
\f[
dS = (r-D)S dt + \sigma S dW
\f]
where \f$ r\f$, \f$ D\f$, and \f$ \sigma \f$ are constants. This is a standard SDE in finance.
The corresponding PDE for the probability density function is 
\f[
\frac{\partial p}{\partial t} =
-\frac{\partial }{\partial S}\{(r-D)S p\}
+ \frac12\frac{\partial^2 }{\partial S^2}\{\sigma^2S^2p\}
\f]
which can also be written as
\f[
\frac{\partial p}{\partial t} =
 (2\sigma^2-r+D)S\frac{\partial p}{\partial S}
+ \frac12\sigma^2S^2\frac{\partial^2 p}{\partial S^2}
+ (\sigma^2-r+D) p.
\f]
So we can reinterpret this in the original form from ParabollicPDE as
\f[
\alpha(S,t) = (2\sigma^2-r+D)S
\f]
\f[
\beta(S,t) = \sigma S
\f]
\f[
\gamma(S,t) = (\sigma^2-r+D)
\f]
For the initial condition, we use 
\f[
f(S) = \begin{cases}
0 \text{ if } S<S_0-1/2\Delta S \text{ or } S>S_0+1/2\Delta S\\
\frac{1}{dS} \text{ if } S_0-1/2\Delta S \leq S \leq S_0+1/2\Delta S
\end{cases} 
\f]
at \f$S=0\f$ we have
\f[
g_1(t) = 0
\f]
and for large \f$S\f$ we use
\f[
g_2(t) = 0.
\f]
Now, given a positive integer \f$n\f$, we can set the parameters up as follows
| Parameter    |   Value      |
|--------------|--------------|
| r            | 0.05         |
| D            | 0.02         |
| sigma        | 0.2          |
| S_0          | 10.0         |
| S_T          | 10.1         |
| theta        | 0.5          |
| t            | 0.0          |
| T            | 1.0          |
| iMax         | n            |
| S_min        | 0.0          |
| S_max        | 5*S_0 = 50   |
| jMax         | 5*n          |
We are able to tun convergence analysis by varying the parameter \f$n\f$. Note that the solution we are looking for is given by the analytic formula for a log normal distribution
\f[
\ln S_T \sim N\left( \ln S_0 + \left(r-D - \frac12 \sigma^2\right)T, \sigma^2 T \right)
\f]
so 
\f[
p(S_T,T) =0.19749617~ .
\f]
Running the numerical solver we get the following results. Note that `p_extrap` is calculated using Richardson extrapolation assuming second order convergence, combining the result from \f$ n \f$ and \f$ n/2 \f$.
|     Steps (n)|      p(S_T,T)|   p_extrap   |      CPU time|
|--------------|--------------|--------------|--------------|
|             8|    0.21027729|   NA         |  (1.6954e-05)|
|            16|    0.19996595|    0.19652884|  (4.1327e-05)|
|            32|    0.19809115|    0.19746622| (0.000152357)|
|            64|    0.19764408|    0.19749506|  (0.00058323)|
|           128|    0.19753311|    0.19749612|  (0.00232197)|
|           256|     0.1975054|    0.19749617|  (0.00926212)|
|           512|    0.19749848|    0.19749617|   (0.0372405)|
|          1024|    0.19749675|    0.19749617|    (0.183431)|

     **/
    int runParabollicPDE_Example1();
    
    
    
    /** @brief Example solving the Black Scholes equation, using ParabollicPDE.
     * @details Solve Black Scholes equation for a European call option.
This is a standard PDE in finance
\f[
\frac{\partial V}{\partial \tau} =
(r-D)S \frac{\partial V}{\partial S}
+ \frac12\sigma^2S^2\frac{\partial^2 V}{\partial S^2} - rV
\f]
where \f$ r\f$, \f$ D\f$, and \f$ \sigma \f$ are constants, and we have used  \f$\tau=T-t\f$ to transform the problem to an initial condition. 

So we can reinterpret this in the original form from ParabollicPDE as
\f[
\alpha(S,\tau) = (r-D)S
\f]
\f[
\beta(S,\tau) = \sigma S
\f]
\f[
\gamma(S,\tau) = -r
\f]
For the initial condition at \f$ \tau=0 \f$, we use 
\f[
f(S) = \begin{cases}
0 \text{ if } S<X\\
S-X \text{ if } S \geq X
\end{cases} 
\f]
at \f$S=0\f$ we have
\f[
g_1(\tau) = 0
\f]
and for large \f$S\f$ we use
\f[
g_2(\tau) = Se^{-D\tau}-Xe^{-r\tau}.
\f]
Now, given a positive integer \f$n\f$, we can set the parameters up as follows
| Parameter    |   Value      |
|--------------|--------------|
| r            | 0.05         |
| D            | 0.02         |
| sigma        | 0.2          |
| S            | 9.375        |
| X            | 10           |
| theta        | 0.5          |
| tau          | 0.0          |
| T            | 1.0          |
| iMax         | n            |
| S_min        | 0.0          |
| S_max        | 5*S_0 = 50   |
| jMax         | 5*n          |
We are able to tun convergence analysis by varying the parameter \f$n\f$. Note that the solution we are looking for is given (in standard notation) as
\f[
V(S,\tau) = Se^{-D\tau}N(d_1)-Xe^{-r\tau}N(d_2)
\f]
so 
\f[
V(S,t=0) = 1.1510188~ .
\f]
Running the numerical solver we get the following results. Note that `V_extrap` is calculated using Richardson extrapolation assuming second order convergence, combining the result from \f$ n \f$ and \f$ n/2 \f$.
|     Steps (n)|      V(S_0)  |   V_extrap   |      CPU time|
|--------------|--------------|--------------|--------------|
|            8 | 1.1267833    | 1.1690444    | (1.4161e-05) |
|           16 | 1.1450066    | 1.1510811    | (3.7572e-05) |
|           32 | 1.1495158    | 1.1510189    | (0.000137554)|
|           64 | 1.1506434    | 1.1510193    | (0.000523982)|
|          128 | 1.150925     | 1.1510188    | (0.002086976)|
|          256 | 1.1509953    | 1.1510188    | (0.008328273)|
|          512 | 1.1510129    | 1.1510188    | (0.033317914)|
|         1024 | 1.1510173    | 1.1510188    | (0.16843009) |
     **/
    int runParabollicPDE_Example2();
    
    /** @brief Example solving Vasicek Bond Option, using ParabollicPDE.
     * @details Consider SDE
     *\f[
dr = \kappa(\theta - r)dt + \sigma dW
\f]
     * then solve for bond, \f$ B(r,t;T_1) \f$
\f[
\frac{\partial B}{\partial \tau} =
\kappa(\theta -r) \frac{\partial B}{\partial r}
+ \frac12\sigma^2\frac{\partial^2 B}{\partial r^2} - rB
\f]
where \f$ \kappa \f$, \f$ \theta\f$, and \f$ \sigma \f$ are constants, and we have used  \f$\tau=T-t\f$ to transform the problem to an initial condition. 

So we can reinterpret this in the original form from ParabollicPDE as
\f[
\alpha(r,\tau) = \kappa(\theta -r)
\f]
\f[
\beta(r,\tau) = \sigma 
\f]
\f[
\gamma(r,\tau) = -r
\f]
For the initial condition at \f$ \tau=0 \f$, we use 
\f[
f(r) = 1 .
\f]
At both boundaries, we assume that the effects of the second order term are minimal, so we solve
\f[
\frac{\partial B}{\partial \tau} =\kappa(\theta -r) \frac{\partial B}{\partial r}-rB.
\f]
Note that this is not the same as saying \f$ \frac{\partial^2 B}{\partial r^2}=0 \f$. So for small values of \f$r\f$ we have the finite difference equation
\f[
\left(\frac{1}{\Delta \tau} +\frac{\kappa(\theta-r_0)}{\Delta r}+0.5r\right)u^{i+1}_0 - \frac{\kappa(\theta-r_0)}{\Delta r}u^{i+1}_1 = \left(\frac{1}{\Delta \tau} -\frac{\kappa(\theta-r_0)}{\Delta r}-0.5r\right)u^{i}_0 + \frac{\kappa(\theta-r_0)}{\Delta r}u^{i}_1
\f]
and for large \f$r\f$ we use
\f[
\left(\frac{1}{\Delta \tau} +\frac{\kappa(\theta-r_\text{jMax})}{\Delta r}+0.5r\right)u^{i+1}_\text{jMax} - \frac{\kappa(\theta-r_\text{jMax})}{\Delta r}u^{i+1}_\text{jMax-1} = \left(\frac{1}{\Delta \tau} -\frac{\kappa(\theta-r_\text{jMax})}{\Delta r}-0.5r\right)u^{i}_\text{jMax} + \frac{\kappa(\theta-r_\text{jMax})}{\Delta r}u^{i}_\text{jMax-1}
\f]
Now, given a positive integer \f$n\f$, we can set the parameters up as follows
| Parameter    |   Value      |
|--------------|--------------|
| r            | 0.05         |
| D            | 0.02         |
| sigma        | 0.2          |
| S            | 9.375        |
| X            | 10           |
| theta        | 0.5          |
| tau          | 0.0          |
| T            | 1.0          |
| iMax         | n            |
| S_min        | 0.0          |
| S_max        | 5*S_0 = 50   |
| jMax         | 5*n          |
We are able to tun convergence analysis by varying the parameter \f$n\f$. Note that the solution we are looking for is given (in standard notation) as
\f[
V(S,\tau) = Se^{-D\tau}N(d_1)-Xe^{-r\tau}N(d_2)
\f]
so 
\f[
V(S,t=0) = 1.1510188~ .
\f]
Running the numerical solver we get the following results. Note that `V_extrap` is calculated using Richardson extrapolation assuming second order convergence, combining the result from \f$ n \f$ and \f$ n/2 \f$.
|     Steps (n)|      V(S_0)  |   V_extrap   |      CPU time|
|--------------|--------------|--------------|--------------|
|            8 | 1.1267833    | 1.1690444    | (1.4161e-05) |
|           16 | 1.1450066    | 1.1510811    | (3.7572e-05) |
|           32 | 1.1495158    | 1.1510189    | (0.000137554)|
|           64 | 1.1506434    | 1.1510193    | (0.000523982)|
|          128 | 1.150925     | 1.1510188    | (0.002086976)|
|          256 | 1.1509953    | 1.1510188    | (0.008328273)|
|          512 | 1.1510129    | 1.1510188    | (0.033317914)|
|         1024 | 1.1510173    | 1.1510188    | (0.16843009) |
     **/
    int runParabollicPDE_Example3();
    
    
    /**  @}   */
    
} 

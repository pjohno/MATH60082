#pragma once
#include <cmath>

namespace MATH60082{
    
    /** @brief Parameters for the Vasicek interest rate model 
     * @details Model is defined as
     *  \f$ dr = \kappa ( \theta - r ) + \sigma_r dW \f$
     * and the variables defined here refer to those parameters
     */
    struct VasicekInterestRateParameters
    {
        /// mean reversion rate
        double kappa;
        /// long term mean
        double theta;
        /// volatility of the process
        double sigma;
    };
    
    /** @brief Functions for the  Vasicek interest rate model 
     * @details Functions allow 
     * 
     * ~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
     * #include "mInterestRateModels.hpp"
     * #include <iostream>
     * using namespace std;
     * using namespace CommonLib::Math;
     * 
     * int main()
     * {
     *   VasicekInterestRateParameters params={0.5,0.05,0.1};
     *   VasicekInterestRateModel V;
     *   double r=0.05,t=0.,T=1.;
     *   // using default values
     *   cout << V.discountBond(r,t,T) << endl;
     *   // set some new parameters
     *   V=params;
     *   cout << V.discountBond(r,t,T) << endl;
     *   // or this way
     *   V={0.25,0.075,0.05};
     *   cout << V.discountBond(r,t,T) << endl;
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~
     * Output from the program is 
     * ~~~~~~~~~~~~~~~~~~~~~~~~~
     * 0.950401
     * 0.952338
     * 0.948823
     * ~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    class VasicekInterestRateModel: protected VasicekInterestRateParameters
    {
        
        double m(double r,double t,double T) const 
        {
            return r*exp(-kappa*(T-t)) + theta*(1-exp(-kappa*(T-t)));
        }
        
        double q(double t,double T) const 
        {
            return sigma*sigma/2./kappa/kappa*(1-exp(-kappa*(T-t)))*(1-exp(-kappa*(T-t)));
        }
        
        double kSquared(double t,double T) const 
        { 
            return sigma*sigma/2./kappa/kappa/kappa*(4.*exp(-kappa*(T-t)) - exp(-2.*kappa*(T-t)) + 2.*kappa*(T-t) - 3.);
        }
        
        double n(double r,double t,double T) const 
        {
            return (T-t)*theta + (r-theta)*(1-exp(-kappa*(T-t)))/kappa;
        }
        
    public:
        /// initialise with some default parameters
        VasicekInterestRateModel();
        /// copy construct from a set of parameters
        VasicekInterestRateModel(const VasicekInterestRateParameters& V);
        /// assign value from a set of parameters
        VasicekInterestRateModel& operator=(const VasicekInterestRateParameters& V);
        /// return the parameters
        const VasicekInterestRateParameters& toVIRP() const;
        
        /** @brief  The mean of dist at time T given r and time t.
         * @details Value is given by 
         * \f[
         * E_t[r_T] = m(r,t,T) - q(t,T)
         * \f]
         * where
         * \f[
         * m(r,t,T) = e^{-\kappa(T-t)}r + (1-e^{-\kappa(T-t)})\theta
         * \f]
         * and
         * \f[
         * q(t,T) = \frac{\sigma^2}{2\kappa^2}(1-e^{-\kappa(T-t)})^2
         * \f]
         */
        double mean(double r,double t,double T) const 
        {
            return m(r,t,T)-q(t,T);
        }
        /** @brief The variance of the distribution at time T given r and time t.
         * @details Value is given by \f[
         * v^2(t,T) = \frac{\sigma^2}{2\kappa}(1-e^{-2\kappa(T-t)})
         * \f]
         */
        double variance(double t,double T) const 
        {
            return (sigma*sigma)/(2.*kappa)*(1-exp(-2.*kappa*(T-t)));
        }
        
        /** @brief Analytic formula for a 1 unit currency discount bond at time t
         * 
         * @details Value is given by \f[
         * \exp\left[\frac12 k^2(t,T) - n(r,t,T)\right]
         * \f]
         * where
         * \f[
         * k^2(t,T) = \frac{\sigma^2}{2\kappa^3}(4e^{-\kappa(T-t)} - e^{-2\kappa(T-t)} + 2\kappa(T-t) - 3)
         * \f]
         * and 
         * \f[
         * n(r,t,T) =(T-t)\theta + (r-\theta)(1-e^{-\kappa(T-t)})/\kappa
         * \f]
         * 
         * @param r is the current interest rate at time t
         * @param t refers to current time t
         * @param T refers to maturity of the bond
         * 
         * @return value of the discount bond
         */
        double discountBond(double r,double t,double T) const 
        {
            return exp(0.5*kSquared(t,T) - n(r,t,T));
        }
        
    };
    
    VasicekInterestRateModel::VasicekInterestRateModel():VasicekInterestRateParameters({0.1,0.1,0.1}){}
    VasicekInterestRateModel::VasicekInterestRateModel(const VasicekInterestRateParameters& V):VasicekInterestRateParameters(V){}
    VasicekInterestRateModel& VasicekInterestRateModel::operator=(const VasicekInterestRateParameters& V){*this = VasicekInterestRateModel(V);return *this;}
    const VasicekInterestRateParameters& VasicekInterestRateModel::toVIRP() const {return *this;}
}

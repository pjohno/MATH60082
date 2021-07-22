
#include <utility>
#include <algorithm>
#include <cmath>

namespace MATH60082
{
    /** \brief bisection algorithm to find the min of a function
     * \details See findMaxBisection for details.
     * @param f the value function, f(x) which is to be minimised.
     * @param xMin minimum value of x, which is lower bound of interval, \f$ a \f$.
     * @param xMax maximum value of x, which is upper bound of interval, \f$ b \f$.
     * @param tol user supplied tolerance, such that \f$ x^m \in [ x^* - tol(b-a) , x^* + tol(b-a) ] \f$ where  \f$ x^m \f$ is the estimate and  \f$ x^* \f$ is the true minimum.
     * @return returns the root and min value
     */
    template<class T,class F>
    static std::pair<T,T> goldenSearch(const F &f,T xMin,T xMax,T tol)
    {
        static const T phi = (3-sqrt(5))/2.;
        // max number of iterations is related to the tolerance through the relation
        const int maxIter = std::max(int(log( tol )/log( 1.-phi )),1);
        // first get value at the left boundary
        // and put in as current min
        T xRoot=xMin; 
        T minValue=f(xRoot);
        //           // next check the +ve golden ratio value
        {
            T x = xMin + phi*(xMax-xMin);
            T temp=f(x);
            if(temp < minValue)
            {
                xRoot=x;
                minValue=temp;
            }
        }
        // the -ve one
        {
            T x = xMax - phi*(xMax-xMin);
            T temp=f(x);
            if(temp < minValue)
            {
                xRoot=x;
                minValue=temp;
            }
        }
        // and finally the right hand side boundary
        {
            T x = xMax;
            T temp=f(x);
            if(temp < minValue)
            {
                xRoot=x;
                minValue=temp;
            }
        }
        // now iterate using golden ratio
        for(int iter=1;iter<=maxIter;iter++)
        {
            // divide the interval [xMin,xMax] into 3 sections using dp
            T h=(xMax-xMin);
            // now evaluate the function at pMin+dp, pMin+2*dp, pMin+3*dp, given that we have already evaluated the function at xMin and xMax previously
            // -- check to see whether anyone is bigger than before
            T x;
            if(xRoot>xMin+0.5*h)
            {
                xMin=xMin+phi*h;
                x = xMax-phi*(xMax-xMin);
            }
            else
            {
                xMax=xMax-phi*h;
                x = xMin+phi*(xMax-xMin);
            }
            T temp=f(x);
            // if this is very close to min value then break
            if(temp < minValue)
            {
                xRoot=x;
                minValue=temp;
            }
            else if(std::fabs(temp-minValue)<tol)
                break;
            
        }
        return {xRoot,minValue};
    }
    
}

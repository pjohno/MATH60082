#pragma once
#include <iostream>
#include <vector>
#include <functional>

namespace MATH60082{
    
    /// @brief a structure to store a finite difference stencil over three points 
    /// @details Assume we have PDE containing the terms
    /// \f[ \alpha(x,t)\frac{\partial v}{\partial x} + \frac12 \beta^2(x,t) \frac{\partial ^2 v}{\partial x^2} \f]
    /// which can be discretised to the form
    /// \f[ a_i~v_{i-1} + b_i~v_i + c_i~v_{i+1} \f]
    /// at the point \f$ (x_\text{min} + i\Delta x , t ) \f$.
    struct Stencil
    {
      /// term multiplying \f[ v_{i-1} \f]
      double a;
      /// term multiplying \f[ v_{i} \f]
      double b;
      /// term multiplying \f[ v_{i+1} \f]
      double c;
    };
    
    /// @brief a structure to store a full finite difference stencil over three points with a right hand side, useful for specifying boundary conditions 
    /// @details Assume we have boundary condition
    /// which can be discretised to the form
    /// \f[ a_i~v_{i-1} + b_i~v_i + c_i~v_{i+1} = d_i\f]
    /// at the point \f$ (x_\text{min} + i\Delta x , t ) \f$.
    struct FullStencil
    {
      /// term multiplying \f[ v_{i-1} \f]
      double a;
      /// term multiplying \f[ v_{i} \f]
      double b;
      /// term multiplying \f[ v_{i+1} \f]
      double c;
      /// term on the right hand side
      double d;
    };
    
    /// return all values in the grid X seperated by a space
    std::ostream& operator<<(std::ostream& output,const Stencil& X);
    /// return all values in the grid X seperated by a space
    std::ostream& operator<<(std::ostream& output,const FullStencil& X);
    
    /// @brief a structure to store a vector of finite difference stencils over an entire grid variable
    /// @details This class contains the functions to create and assign a set of finite difference stencils given 
    /// an SDE. We do allow manual updating of the stencils, but not resizing.
    class VectorStencil
    {
    protected:
      typedef unsigned int size_t;
      /// storage for the stencils
      std::vector<Stencil> S;
      
    public:
      
      /// construct a set of empty stencils
      VectorStencil(size_t n=0):S(n){}
      
      /// create stencils for the grid x at time t
      /// @details Assuming we have PDE containing the terms
      /// \f[ \alpha(x,t)\frac{\partial v}{\partial x} + \frac12 \beta^2(x,t) \frac{\partial ^2 v}{\partial x^2} \f]
      /// we choose either central or one-sided derivatives for the first order terms to make sure that \f$ a > 0\f$, \f$ b <0 \f$ and \f$ c >0\f$
      /// to maintain stability.
      virtual int createStencils(const std::function<double(double,double)> &alpha,const std::function<double(double,double)> &beta,const std::vector<double>& x,const double &t=0.);
      /// reset stencil
      void resetStencil(size_t i,const Stencil& s=Stencil()){S[i]=s;}
      
      /// return an read write reference to the stencil at index i
      Stencil& operator[](size_t index){return S[index];}
      /// return a read only reference to the stencil at index i
      const Stencil& operator[](size_t index) const {return S[index];}
      /// return the size of the vector
      size_t size(void) const {return S.size();}
      
    };
    
    /// return all values in the grid X seperated by a space
    std::ostream& operator<<(std::ostream& output,const VectorStencil& X);
    
}

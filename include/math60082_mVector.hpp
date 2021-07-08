#pragma once

#include <iostream>
// #define _MATH60082_DEBUG__

namespace MATH60082{
    
    /** @brief MVector is a mathematical vector that can has the normal mathematical operations available
     * @details This will interface with matrices to provide 
     */
    class MVector
    {
        unsigned int N;
        double* v;
    public:
        explicit MVector():N(0),v(nullptr){/*std::cout << "intialise empty\n";*/}
        explicit MVector(unsigned int n):N(n),v(new double[n]){
            #ifdef _MATH60082_DEBUG__ 
            std::cout << "intialise NEW " << v << " \n";
            #endif
        }
        explicit MVector(unsigned int n,double x);
        
        ~MVector(){
            #ifdef _MATH60082_DEBUG__ 
            std::cout << "delete " << v << "\n";
            #endif
            if(v!=nullptr) delete [] v;
            
        }
        // normal copies
        MVector(const MVector& X):N(X.N),v(new double[X.N]){std::copy(X.v,X.v+X.N,v);
            #ifdef _MATH60082_DEBUG__ 
            std::cout << " copy construct NEW " << v << " " << X.v << "\n";
            #endif
        }
        MVector& operator=(const MVector &X);
        // c++11 initialiser and move construct, move=
        explicit MVector(std::initializer_list<double> list);
        MVector(MVector&& X):N(X.N),v(X.v){
            #ifdef _MATH60082_DEBUG__ 
            std::cout << "move construct " << v<< " " << X.v  << "\n";
            #endif
            X.N=0;X.v=nullptr;
        }
        MVector& operator=(MVector&& X){
            if(v!=nullptr){
                if(v==X.v)
                {
                    #ifdef _MATH60082_DEBUG__ 
                    std::cout << " do nothing " << v << " " << X.v << "\n";
                    #endif	    
                    
                    return *this;
                }
                #ifdef _MATH60082_DEBUG__ 
                std::cout << "delete " << v << " ";
                #endif
                delete [] v;
                
            }
            #ifdef _MATH60082_DEBUG__ 
            std::cout << v << "= &&";
            #endif
            N=X.N;v=X.v;X.N=0;X.v=nullptr;
            #ifdef _MATH60082_DEBUG__ 
            std::cout << v << "\n";
            #endif
            return *this;
        }
        
        double& operator[](int index){return v[index];}
        double operator[](int index) const {return v[index];}
        unsigned int size() const {return N;}
        // operations on a vector
        double maxNorm() const;
        double l2Norm() const;
        // resize vector
        void resize(int n);
        void resize(int n,double x);
        // return array
        double* returnArray(){return v;}
        const double* returnArray() const {return v;}
    };
    
    std::ostream& operator<<(std::ostream &output,const MVector &X);
    
    // normal scalar mult
    MVector operator*(const double& lhs,const MVector &rhs);
    MVector operator*(const MVector& lhs,const double &rhs);
    // c++11 scalar mult
    MVector operator*(const double& lhs,MVector&& rhs);
    // normal addition
    MVector operator+(const MVector& lhs,const MVector &rhs);
    // c++11 additions
    MVector operator+(const MVector& lhs,MVector&& rhs);
    MVector operator+(MVector&& lhs,const MVector& rhs);
    // normal subtraction
    MVector operator-(const MVector& lhs,const MVector &rhs);
    // c++11 subtractions
    MVector operator-(const MVector& lhs,MVector&& rhs);
    MVector operator-(MVector&& lhs,const MVector &rhs);
    
    MVector operator/(const MVector& lhs,const double &rhs);
    
    double dot(const MVector& lhs,const MVector& rhs);
    
    
}


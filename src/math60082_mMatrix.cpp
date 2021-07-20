#include "math60082_mMatrix.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace MATH60082
{
    
    MMatrix::MMatrix():N(0),M(0),A(nullptr){
        #ifdef DEBUG_MMATRIX_CONSTRUCT
        std::cout << " Construct empty MMatrix(int "<<N<<",int "<<M<<") " << this <<" -- A="<<A<<std::endl;
        #endif 
        
    }
    
    MMatrix::MMatrix(const MMatrix &X):N(X.N),M(X.M),A(new double [X.N*X.M])
    {
        std::copy(X.A,X.A+X.size(),A);
        #ifdef DEBUG_MMATRIX_CONSTRUCT
        std::cout << " Copy Construct MMatrix(const MMatrix "<<&X<<") " << this <<" -- A="<<A<<std::endl;
        #endif 
        
    }
    
    MMatrix::MMatrix(std::initializer_list<std::initializer_list<double> > list):N(list.size()),M(0),A(nullptr)
    {
        for(auto &row : list)
        {
            M = std::max((std::size_t)M,row.size());
        }
        A = new double [N*M];
        unsigned int i=0;
        for(auto &row : list)
        {
            unsigned int j=0;
            for(auto &element : row)
            {
                if(j>=M)throw;
                A[i*M+j]=element;
                j++;
            }
            for(unsigned int jStar=j;jStar<M;jStar++)
            {
                A[i*M+j]=0.;
            }
            i++;
        }
        #ifdef DEBUG_MMATRIX_CONSTRUCT
        std::cout << " Copy Construct MMatrix(std::initializer_list<std::initializer_list<double> > list) " << this <<" -- A="<<A<<std::endl;
        #endif 
        
    }
    
    MVector MMatrix::operator*(const MVector &x) const 
    {
        MVector temp(rows());
        double *posA=A;
        double *tempA=temp.returnArray();
        for(unsigned int i=0;i<N;i++)
        {
            double sum=0.;
            const double *posx=x.returnArray();
            for(unsigned int j=0;j<M;j++)sum+=*posA++ * *posx++;
            tempA[i] = sum;
        }
        return temp;
    }
    
    MMatrix& MMatrix::operator=(const MMatrix &X)
    {
        if(&X==this)
            return *this;
        resize(X.N,X.M);
        std::copy(X.A,X.A+X.size(),A);
        #ifdef DEBUG_MMATRIX
        std::cout << " Copy = MMatrix(int "<<N<<",int "<<M<<"):A="<<A<<" <=="<<X.A<<std::endl;
        #endif 
        
        return *this;
    }
    
    MMatrix& MMatrix::operator=(const double &x)
    {
        for(unsigned int i=0;i<size();i++)A[i]=x;
        return *this;
    }
    
    std::ostream& operator<<(std::ostream &output,const MMatrix &X)
    {
        
        unsigned int imorethanone=std::min((unsigned int)1,X.rows());
        for(unsigned int i=0;i<imorethanone;i++)
        {
            output << "( ";
            unsigned int jmorethanone=std::min((unsigned int)1,X.cols());
            for(unsigned int j=0;j<jmorethanone;j++)
                output << X[i][j];
            for(unsigned int j=jmorethanone;j<X.cols();j++)
                output << " " << X[i][j];
            output << " )\n";
        }
        for(unsigned int i=imorethanone;i<X.rows();i++)
        {
            output << "( ";
            unsigned int jmorethanone=std::min((unsigned int)1,X.cols());
            for(unsigned int j=0;j<jmorethanone;j++)
                output << X[i][j];
            for(unsigned int j=jmorethanone;j<X.cols();j++)
                output << " " << X[i][j];
            output << " )\n";
        }
        return output;
    }
    
    void MMatrix::resize(unsigned int n,unsigned int m)
    {
        #ifdef DEBUG_MMATRIX
        std::cout << " resize MMatrix("<<N<<","<<M<<")==>("<<n<<","<<m<<")\n";
        #endif 
        
        if(N!=n || M!=m)
            clear();
        else
            return;
        A = new double[n*m];
        N=n;M=m;
    }
    
    void MMatrix::resize(unsigned int n,unsigned int m,double x)
    {
        resize(n,m);
        for(unsigned int i=0;i<size();i++)
            A[i]=x;
    }
    
    void MMatrix::clear()
    {
        N=0;
        M=0;
        if(A!=nullptr)
            delete [] A;
        A=nullptr;
    }
    
    double MMatrix::l2Norm() const
    {
        double temp=0.;
        const double *tempPos=returnArray();
        for(unsigned int i=0;i<size();i++)
        {
            temp+=(*tempPos)*(*tempPos);
            tempPos++;
        }
        return sqrt(temp);
    }
    
    MMatrix operator*(const double& lhs,const MMatrix &rhs)
    {
        MMatrix temp(rhs);
        double *tempPos=temp.returnArray();
        for(unsigned int i=0;i<rhs.size();i++)
            (*tempPos++)*=lhs;
        return temp;
    }
    
    MMatrix operator*(const double& lhs,MMatrix&& rhs)
    {
        MMatrix temp(std::move(rhs));
        double *tempPos=temp.returnArray();
        for(unsigned int i=0;i<temp.size();i++)
            (*tempPos++)*=lhs;
        return temp;
    }
    
    MMatrix operator+(const MMatrix& lhs,const MMatrix &rhs)
    {
        MMatrix temp(lhs);
        double *tempPos=temp.returnArray();
        const double *rPos=rhs.returnArray();
        for(unsigned int i=0;i<lhs.size();i++)
            (*tempPos++)+=(*rPos++);
        return temp;
    }
    
    MMatrix operator+(const MMatrix& lhs,MMatrix&& rhs)
    {
        MMatrix temp(std::move(rhs));
        double *tempPos=temp.returnArray();
        const double *lPos=lhs.returnArray();
        for(unsigned int i=0;i<lhs.size();i++)
            (*tempPos++)+=(*lPos++);
        return temp;
    }
    
    MMatrix operator+(MMatrix&& lhs,const MMatrix& rhs)
    {
        MMatrix temp(std::move(lhs));
        double *tempPos=temp.returnArray();
        const double *rPos=rhs.returnArray();
        for(unsigned int i=0;i<rhs.size();i++)
            (*tempPos++) +=(*rPos++);
        return temp;
    }
    
    MMatrix operator-(const MMatrix& lhs,const MMatrix &rhs)
    {
        MMatrix temp(lhs);
        double *tempPos=temp.returnArray();
        const double *rPos=rhs.returnArray();
        for(unsigned int i=0;i<lhs.size();i++)
            (*tempPos++)-=(*rPos++);
        return temp;
    }
    
    
    MMatrix operator-(const MMatrix& lhs,MMatrix&& rhs)
    {
        MMatrix temp(std::move(rhs));
        double *tempPos=temp.returnArray();
        const double *lPos=lhs.returnArray();
        for(unsigned int i=0;i<lhs.size();i++)
        {
            *tempPos*=-1;
            (*tempPos++)+=(*lPos++);
        }
        return temp;
    }
    
    MMatrix operator-(MMatrix&& lhs,const MMatrix& rhs)
    {
        MMatrix temp(std::move(lhs));
        double *tempPos=temp.returnArray();
        const double *rPos=rhs.returnArray();
        for(unsigned int i=0;i<rhs.size();i++)
            (*tempPos++)-=(*rPos++);
        return temp;
    }
    
    double dot(const MMatrix& lhs,const MMatrix& rhs)
    {
        double temp=0.;
        const double *lPos=lhs.returnArray();
        const double *rPos=rhs.returnArray();
        for(unsigned int i=0;i<lhs.size();i++)
            temp+=(*rPos++)*(*lPos++);
        return temp;
    }
    
}

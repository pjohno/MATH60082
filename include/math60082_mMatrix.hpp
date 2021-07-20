#pragma once

#include <iostream>
#include <ostream>
#include "math60082_mVector.hpp"

namespace MATH60082 
{
    
    /** @brief MMatrix contains an array of doubles, that can be accessed as if they were a matrix with (i,j) operator
     * @details The MMatrix has access operators and mathematical operations (add, subtract, etc.) that make it more like a matrix in matlab.
     */
    class MMatrix
    {
        /// storage for the number of rows in the matrix
        unsigned int N;
        /// storage for the number of columns in the matrix
        unsigned int M;
        /// Pointer to the storage for the array. We use pointers here to allow std::move constructs minimizing unnecessary copies
        double* A;
        
    public:
        /// default constructor
        explicit MMatrix();
        /// construct a matrix with n rows and m columns
        explicit MMatrix(int n,int m):N(n),M(m),A(new double [n*m])
        {
            #ifdef DEBUG_MMATRIX_CONSTRUCT
            std::cout << " Construct MMatrix(int "<<n<<",int "<<m<<") " << this <<" -- A="<<A<<std::endl;
            #endif 
        }
        /// construct a matrix with n rows and m columns and assign x to all elements in the matrix
        explicit MMatrix(int n,int m,double x):N(n),M(m),A(new double [n*m])
        {
            for(unsigned int i=0;i<size();i++)A[i]=x;
            #ifdef DEBUG_MMATRIX_CONSTRUCT
            std::cout << " Construct MMatrix(int "<<n<<",int "<<m<<",double "<<x<<") " << this <<" -- A="<<A<<std::endl;
            #endif 
        }
        /// construct the matrix making a new version and copying the data 
        MMatrix(const MMatrix &X);
        /// construct a new matrix by moving the data, here matrix X will be deleted and no longer usable
        explicit MMatrix(MMatrix&& X):N(X.N),M(X.M),A(X.A)
        {
            X.N=0;X.M=0;X.A=nullptr;
            #ifdef DEBUG_MMATRIX_CONSTRUCT
            std::cout << " Copy Contruct MMatrix(MMatrix&& "<<&X<<") " << this <<" -- A="<<A<<std::endl;;
            #endif
        }
        /** @brief construct from a double list
         * @details An example here would be 
         * ~~~{.c}
         * MMatrix A_list_initializer({ { 1. , 2. , 3.} , { 4. , 5. , 6.}  , { 7. , 8. , 9.} });
         * cout << A_list_initializer;
         * ~~~
         * and output would be 
         * ~~~
         * ( 1 2 3 )
         * ( 4 5 6 )
         * ( 7 8 9 )
         * ~~~
         * The code 
         * ~~~{.c}
         * MMatrix A_columns_different ({ { 1. , 2. , 3. } , { 4. , 5. } , { 6. , 7. } , { 8. , 9. , 10. , 11. } });
         * cout << A_columns_different;
         * ~~~
         * will output
         * ~~~
         * ( 1 2 3 0 )
         * ( 4 5 0 0 )
         * ( 6 7 0 0 )
         * ( 8 9 10 11 )
         * ~~~
         */
        explicit MMatrix(std::initializer_list<std::initializer_list<double> > list);
        /// will assign the matrix in current storage (if conformal, deletes/creates new storage if not), copying the data over 
        MMatrix& operator=(const MMatrix &X);
        /// will update all elements in the matrix \f$ a_{ij}=x \f$
        MMatrix& operator=(const double &x);
        /** @brief will assign the current matrix to the temporary one created
         * @details This operator allows additions/subtraction/multiplications to chain together by only ever creating one temporary matrix.
         * 
         * The code
         * ~~~{.c}
         * MMatrix A( { { 1. , 2. } , { 3. , 4.} } ),B( { { 5. , 6. } , { 7. , 8.} } );
         * cout << "A=\n"<< A << "B=\n" << B << endl;
         * MMatrix C(A);
         * C(1,1) = 1.;
         * cout << "C=\n"<< C << endl;
         * A = A + (B + 3*C);
         * cout << "A=\n"<< A << "B=\n" << B << "C=\n" << C << endl;
         * ~~~
         * only ever create one termporary variable. Define DEBUG_MMATRIX_CONSTRUCT to see this in action. The output would look like this below:
         * ~~~
         *  Copy Construct MMatrix(std::initializer_list<std::initializer_list<double> > list) 0x7ffeacca2060 -- A=0x11c90f0
         * Copy Construct MMatrix(std::initializer_list<std::initializer_list<double> > list) 0x7ffeacca2070 -- A=0x11c9120
         * A=
         * ( 1 2 )
         * ( 3 4 )
         * B=
         * ( 5 6 )
         * ( 7 8 )
         * 
         * Copy Construct MMatrix(const MMatrix 0x7ffeacca2060) 0x7ffeacca2080 -- A=0x11c9150
         * C=
         * ( 1 2 )
         * ( 3 1 )
         * 
         * Copy Construct MMatrix(const MMatrix 0x7ffeacca2080) 0x7ffeacca2090 -- A=0x11c9180
         * Copy Contruct MMatrix(MMatrix&& 0x7ffeacca2090) 0x7ffeacca20a0 -- A=0x11c9180
         * Copy Contruct MMatrix(MMatrix&& 0x7ffeacca20a0) 0x7ffeacca20b0 -- A=0x11c9180
         * MMatrix& operator=(MMatrix&& X) -- delete A=0x11c90f0
         * MMatrix& operator=(MMatrix&& X) -- swap A=0x11c9180
         * delete MMatrix 0x7ffeacca20b0 -- delete A=0
         * delete MMatrix 0x7ffeacca20a0 -- delete A=0
         * delete MMatrix 0x7ffeacca2090 -- delete A=0
         * A=
         * ( 9 14 )
         * ( 19 15 )
         * B=
         * ( 5 6 )
         * ( 7 8 )
         * C=
         * ( 1 2 )
         * ( 3 1 )
         * 
         * delete MMatrix 0x7ffeacca2080 -- delete A=0x11c9150
         * delete MMatrix 0x7ffeacca2070 -- delete A=0x11c9120
         * delete MMatrix 0x7ffeacca2060 -- delete A=0x11c9180
         * 
         * ~~~
         */
        MMatrix& operator=(MMatrix&& X)
        {
            if(N==X.N && M==X.M && A==X.A)
                return *this;
            else
            {
                #ifdef DEBUG_MMATRIX_CONSTRUCT
                std::cout << " MMatrix& operator=(MMatrix&& X) -- delete A="<<A<<std::endl;;
                #endif
                delete [] A;
            }
            N=X.N;M=X.M;A=X.A;
            X.N=0;X.M=0;X.A=nullptr;
            #ifdef DEBUG_MMATRIX_CONSTRUCT
            std::cout << " MMatrix& operator=(MMatrix&& X) -- swap A="<<A<<std::endl;
            #endif
            return *this;
        }
        
        /** @brief default deconstructor
         */
        virtual ~MMatrix(){
            #ifdef DEBUG_MMATRIX_CONSTRUCT
            std::cout << " delete MMatrix "<<this<<" -- delete A="<<A<<std::endl;
            #endif
            if(A!=nullptr) delete [] A;
        }
        
        /** @brief column access
         *   @details allow read only access A[i][j] for MMatrix A
         *   @return pointer to the first element in the \f$ i \f$ th row. 
         */
        const double* operator[](unsigned int i) const {return &(A[i*M]);}
        /** @brief column access
         *   @details allow read/write access A[i][j] for MMatrix A
         *   @return pointer to the first element in the \f$ i \f$ th row. 
         */
        double* operator[](unsigned int i){return &(A[i*M]);}
        //       /** @brief element access
        //        *   @details allow read only access A(i,j) for MMatrix A
        //        *   @return value of the element in the \f$ i \f$ th row and \f$ j \f$ th column. 
        //        */
        //       double operator()(unsigned int i,unsigned int j) const {return A[i*M+j];}
        //       /** @brief element access
        //        *   @details allow read/write access A(i,j) for MMatrix A
        //        *   @return value of the element in the \f$ i \f$ th row and \f$ j \f$ th column. 
        //        */
        //       double& operator()(unsigned int i,unsigned int j){return A[i*M+j];}
        /** @brief size access
         *   @details return the number of rows
         *   @return number of rows \f$ N \geq 0 \f$
         */
        unsigned int rows() const {return N;}
        /** @brief size access
         *   @details return the number of columns
         *   @return number of columns \f$ M \geq 0 \f$
         */
        unsigned int cols() const {return M;}
        /** @brief size access
         *   @details total number of elements, also the size of the array storage
         *   @return number of total number of elements \f$ N*M \geq 0 \f$
         */
        unsigned int size() const {return N*M;}
        /** @brief Resize matrix
         *   @details allow the matrix to be resized
         *   @param n \f$ n \geq 0 \f$ is the number of rows in the new matrix
         *   @param m \f$ m \geq 0 \f$ is the number of columns in the new matrix
         */
        void resize(unsigned int n,unsigned int m);
        /** @brief Resize matrix
         *   @details allow the matrix to be resized
         *   @param n \f$ n \geq 0 \f$ is the number of rows in the new matrix
         *   @param m \f$ m \geq 0 \f$ is the number of columns in the new matrix
         *   @param x is the default value for each element in the new matrix \f$ a_{ij}=x \f$ 
         */
        void resize(unsigned int n,unsigned int m,double x);
        /// clear the matrix
        void clear();
        /** @brief Return the l2 norm of the matrix
         *   @details This is calcuted by the following
         *   \f[
         *    L = \sqrt{ \sum_i \sum_j a_{ij}^2 }
         *   \f]
         * @return L as above
         */
        double l2Norm() const;
        
        /// Return a read/write access pointer to the first element in the array storage
        double* returnArray(){return A;}
        /// Return pointer to the array storage
        double** returnArrayPtr(){return &A;}
        /// Return a read only pointer to the first element in the array storage
        const double* returnArray() const {return A;}
        
        /// Implement a matrix/vector multiplication 
        MVector operator*(const MVector &x) const ;
        
        /// Normal scalar mult
        /// This will copy/construct a new matrix in which to put the result
        friend MMatrix operator*(const double& lhs,const MMatrix &rhs);
        /// Move scalar mult
        /// This will move the result into the matrix rhs, no new matrices created
        friend MMatrix operator*(const double& lhs,MMatrix&& rhs);
        
        
        /** @brief  Normal addition
         * @details  This will copy/construct a new matrix in which to put the result.
         * @param A a matrix \f$ A \f$
         * @param B a matrix \f$ B \f$ with the same rows/columns as \f$ A \f$
         * @return \f$ A+B \f$
         */
        friend MMatrix operator+(const MMatrix& A,const MMatrix &B);
        /** @brief  Move additions
         * @details  This will move the result into the matrix B, no new matrices created
         * @param A a matrix \f$ A \f$
         * @param B a matrix \f$ B \f$ with the same rows/columns as \f$ A \f$, must be a temporary matrix not used again
         * @return \f$ A+B \f$
         */
        friend MMatrix operator+(const MMatrix& A,MMatrix&& B);
        /** @brief  Move additions
         * @details This will move the result into the matrix A, no new matrices created
         * @param A a matrix \f$ A \f$, must be a temporary matrix not used again
         * @param B a matrix \f$ B \f$ with the same rows/columns as \f$ A \f$
         * @return \f$ A+B \f$
         */
        friend MMatrix operator+(MMatrix&& A,const MMatrix& B);
        /** @brief normal subtraction
         * @details This will copy/construct a new matrix in which to put the result
         * @param A a matrix \f$ A \f$
         * @param B a matrix \f$ B \f$ with the same rows/columns as \f$ A \f$
         * @return \f$ A-B \f$
         */
        friend MMatrix operator-(const MMatrix& A,const MMatrix &B);
        /** @brief  Move subtractions
         * @details This will move the result into the matrix B, no new matrices created
         * @param A a matrix \f$ A \f$
         * @param B a matrix \f$ B \f$ with the same rows/columns as \f$ A \f$, must be a temporary matrix not used again
         * @return \f$ A-B \f$
         */
        friend MMatrix operator-(const MMatrix& A,MMatrix&& B);
        /** @brief  Move subtractions
         * @details This will move the result into the matrix A, no new matrices created
         * @param A a matrix \f$ A \f$, must be a temporary matrix not used again
         * @param B a matrix \f$ B \f$ with the same rows/columns as \f$ A \f$
         * @return \f$ A-B \f$
         */
        friend MMatrix operator-(MMatrix&& A,const MMatrix &B);
        
        /** @brief calculate the dot product of two MMatrices as if they were rolled out vectors
         * @details This dot product assumes that the matrices have been rolled out into vectors so that
         * \f[ v_{iM+j}=a_{ij}\f]
         * and
         * \f[ \bm{A.B} = \sum_i \sum_j a_{ij} b_{ij} \f]
         * @return \f$ \bm{A}.\bm{B} \f$
         */
        friend double dot(const MMatrix& A,const MMatrix& B);
        
        /** @brief Output the matrix in a standard format
         * @details Output an (n+1) x (m+1) matrix in format
         * ~~~
         * ( a_00 a_01 a_02 ... a_0m ) 
         * ( a_10 a_11 a_12 ... a_1m ) 
         * ( .... .... .... ... .... ) 
         * ( a_n0 a_n1 a_n2 ... a_nm ) 
         * ~~~
         * @param X matrix to be printed
         */
        friend std::ostream& operator<<(std::ostream& output,const MMatrix& X);
        
    }; // end class MMatrix
    
}

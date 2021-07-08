#include "math60082_mVector.hpp"
#include <cmath>
#include <algorithm>

namespace MATH60082 {

std::ostream& operator<<(std::ostream &output,const MVector &X)
{
  //output << "( ";
  unsigned int morethanone=std::min((unsigned int)(1),X.size());
  for(unsigned int i=0;i<morethanone;i++)
    output << X[i];
  for(unsigned int i=morethanone;i<X.size();i++)
    output << " , " << X[i];
  //output << " )";
  return output;
}

MVector::MVector(std::initializer_list<double> list):N(list.size()),v(new double[list.size()])
{
  #ifdef _MY_DUBUG__ 
  std::cout << "intialise NEW " << N << " " << v << "\n";
  #endif
  auto  vPos=v;for(auto lPos=list.begin();lPos!=list.end();lPos++,vPos++)*vPos=*lPos;/*std::cout << "intialise list\n";*/
}

MVector::MVector(unsigned int n,double x):N(n),v(new double[n])
{
  for(unsigned int i=0;i<N;i++)v[i]=x;
  #ifdef _MY_DUBUG__ 
  std::cout << "intialise NEW " << n << " " << v << "\n";
  #endif
}

void MVector::resize(int n,double x)
{
  resize(n);
  double *vPos=v;
  for(unsigned int i=0;i<N;i++)*vPos++=x;
}

void MVector::resize(int n)
{
  #ifdef _MY_DUBUG__ 
  std::cout << "delete " << v << "\n";
  #endif
  if(v!=nullptr) delete [] v;
  
  v=new double[n];
  N=n;
  #ifdef _MY_DUBUG__ 
  std::cout << "resize NEW " << n << " " << v << "\n";
  #endif
}

double MVector::maxNorm() const
{
  double temp(0.);
  if(size()>0)temp=std::abs(v[0]);
  for(unsigned int i=1;i<size();i++)temp=std::max(temp,v[i]);
  return temp;
}
double MVector::l2Norm() const
{
  double temp=0.;
  for(unsigned int i=0;i<size();i++)temp+=v[i]*v[i];
  return sqrt(temp);
}

MVector operator*(const double& lhs,const MVector &rhs)
{
  MVector temp(rhs);
  for(unsigned int i=0;i<rhs.size();i++)
    temp[i]*=lhs;
  return temp;
}

MVector operator*(const double& lhs,MVector&& rhs)
{
  MVector temp(std::move(rhs));
  double *tempPos=temp.returnArray();
  for(unsigned int i=0;i<temp.size();i++)
    (*tempPos++)*=lhs;
  return temp;
}

MVector operator*(const MVector& lhs,const double &rhs)
{
  MVector temp(lhs);
  for(unsigned int i=0;i<lhs.size();i++)
    temp[i]*=rhs;
  return temp;
}

MVector operator/(const MVector& lhs,const double &rhs)
{
  MVector temp(lhs);
  for(unsigned int i=0;i<lhs.size();i++)
    temp[i]/=rhs;
  return temp;
}

MVector operator+(const MVector& lhs,const MVector &rhs)
{
  MVector temp(lhs);
  double *tempPos=temp.returnArray();
  const double *rPos=rhs.returnArray();
  for(unsigned int i=0;i<lhs.size();i++)
    (*tempPos++)+=(*rPos++);
  return temp;
}

MVector operator+(const MVector& lhs,MVector&& rhs)
{
  MVector temp(std::move(rhs));
  double *tempPos=temp.returnArray();
  const double *lPos=lhs.returnArray();
  for(unsigned int i=0;i<lhs.size();i++)
    (*tempPos++)+=(*lPos++);
  return temp;
}

MVector& MVector::operator=(const MVector &X)
{
  #ifdef _MY_DUBUG__ 
  std::cout << v<< " ";
  #endif
  if(v!=nullptr && N!=X.N)
  {
    #ifdef _MY_DUBUG__ 
    std::cout << " delete " << v << " ";
    #endif
    delete [] v;
  }
  if(N!=X.N)
  {
    #ifdef _MY_DUBUG__ 
    std::cout << " NEW ";
    #endif
    v = new double [X.N];
  }
  #ifdef _MY_DUBUG__ 
  std::cout << v << " = & " << X.v << "\n";
  #endif
  N = X.N;
  std::copy(X.v,X.v+X.N,v);
  return *this;
}

MVector operator-(const MVector& lhs,const MVector &rhs)
{
  MVector temp(lhs);
  for(unsigned int i=0;i<lhs.size();i++)
    temp[i]-=rhs[i];
  return temp;
}

MVector operator-(const MVector& lhs,MVector&& rhs)
{
  MVector temp(std::move(rhs));
  double *tempPos=temp.returnArray();
  const double *lPos=lhs.returnArray();
  for(unsigned int i=0;i<lhs.size();i++)
  {
    *tempPos*=-1;
    (*tempPos++)+=(*lPos++);
  }
  return temp;
}

double dot(const MVector& lhs,const MVector& rhs)
{
  double temp=0.;
  for(unsigned int i=0;i<lhs.size();i++)
    temp+=rhs[i]*lhs[i];
  return temp;
}

}



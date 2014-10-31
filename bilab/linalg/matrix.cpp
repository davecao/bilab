#include "matrix.h"

using namespace std;

template <class T>
inline void CWmatrix<T>::show() const
{
   int width = 10;
   string str;
   cout.precision(4);

   cout<<endl;
   cout<<str.insert(0, mm, '-')<<endl;
  for (int i=0; i<nn; i++) {
    for (int j=0; j<mm; j++) {
      cout<<setw(width)<<v[i][j];
    }
    cout<<endl;
  }
  cout<<endl;
}

template <class T>
inline void CWmatrix<T>::show(string info) const
{
   int width = 10;
   int l = info.length();
   string str;
   cout.precision(4);

   cout<<endl;
   cout<<info<<endl;
   cout<<str.insert(0, l, '-')<<endl;
  for (int i=0; i<nn; i++) {
    for (int j=0; j<mm; j++) {
      cout<<setw(width)<<v[i][j];
    }
    cout<<endl;
  }
  cout<<endl;
}

template <class T>
CWmatrix<T>::CWmatrix() : nn(0), mm(0), v(NULL) {}

template <class T>
CWmatrix<T>::CWmatrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL){
  int i,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
CWmatrix<T>::CWmatrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
  int i,j,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< n; i++) v[i] = v[i-1] + m;
  for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
CWmatrix<T>::CWmatrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
   int i,j,nel=m*n;
   if (v) v[0] = nel>0 ? new T[nel] : NULL;
   for (i=1; i< n; i++) v[i] = v[i-1] + m;
   for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
CWmatrix<T>::CWmatrix(const CWmatrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
  int i,j,nel=mm*nn;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
CWmatrix<T> & CWmatrix<T>::operator=(const CWmatrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//    if matrix and rhs were different sizes, matrix
//    has been resized to match the size of rhs
{
  if (this != &rhs) {
    int i,j,nel;
    if (nn != rhs.nn || mm != rhs.mm) {
      if (v != NULL) {
        delete[] (v[0]);
        delete[] (v);
      }
      nn=rhs.nn;
      mm=rhs.mm;
      v = nn>0 ? new T*[nn] : NULL;
      nel = mm*nn;
      if (v) v[0] = nel>0 ? new T[nel] : NULL;
      for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
    }
    for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
  }
  return *this;
}
template<class T>
CWmatrix<T> & CWmatrix<T>::operator=(const T& a){
   //Obtain array size
   for (int i=0; i< nn; i++) for (int j=0; j<mm; j++) v[i][j] = a;
}
template <class T>
inline T* CWmatrix<T>::operator[](const int i)  {//subscripting: pointer to row i
#ifdef _CHECKBOUNDS_
   if (i<0 || i>=nn) {
      throw("Matrix subscript out of bounds");
   }
#endif
   return v[i];
}

template <class T>
inline const T* CWmatrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
   if (i<0 || i>=nn) {
      throw("Matrix subscript out of bounds");
   }
#endif
   return v[i];
}

template <class T>
inline T& CWmatrix<T>::operator()(const int i,const int j)  {//subscripting: pointer to row i
#ifdef _CHECKBOUNDS_
  if (i<0 || i>=nn) {
    throw("Matrix: Row's subscript out of bounds");
  }
  if (j<0 || j>=mm) {
    throw("Matrix: Column's subscript out of bounds");
  }
#endif
  return v[i][j];
}

template <class T>
inline const T& CWmatrix<T>::operator()(const int i,const int j) const
{
#ifdef _CHECKBOUNDS_
  if (i<0 || i>=nn) {
    throw("CWmatrix: Row's subscript out of bounds");
  }
  if (j<0 || j>=mm) {
    throw("CWmatrix: Column's subscript out of bounds");
  }
#endif
  return v[i][j];
}


template <class T>
inline int CWmatrix<T>::nrows() const
{
  return nn;
}

template <class T>
inline int CWmatrix<T>::ncols() const
{
  return mm;
}

template <class T>
void CWmatrix<T>::resize(int newn, int newm)
{
  int i,nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new T*[nn] : NULL;
    nel = mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  }
}

template <class T>
void CWmatrix<T>::assign(int newn, int newm, const T& a)
{
  int i,j,nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new T*[nn] : NULL;
    nel = mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  }
  for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}

template <class T>
CWmatrix<T> operator+(const CWmatrix<T> a, const CWmatrix<T> b)
{
  int rows = a.nrows();
  int cols = a.ncols();

  CWmatrix<T> temp(rows,cols);
  if((a.nrows() == b.nrows()) && (a.ncols() == b.ncols()))
  {
    for(int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
        temp[i][j] = a[i][j] + b[i][j];
      }
    }
  }else{
    throw("Size of Matrice-> not same");
  }
  return temp;
}

template <class T>
CWmatrix<T> operator-(const CWmatrix<T> a, const CWmatrix<T> b)
{
  int rows = a.nrows();
  int cols = a.ncols();

  CWmatrix<T> temp(rows,cols);
  if((a.nrows() == b.nrows()) && (a.ncols() == b.ncols()))
  {
    for(int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
        temp[i][j] = a[i][j] - b[i][j];
      }
    }
  }else{
    throw("Matrice are not in the same size.");
  }
  return temp;
}

template <class T>
CWmatrix<T> operator*(const CWmatrix<T> a, const CWmatrix<T> b)
{
   int i,j,k,nthreads,chunk;
   int nra = a.nrows();
   int nca = a.ncols();
   int ncb = b.ncols();
   int nrb = b.nrows();

   if(nca != nrb){
      throw("Matrix size INCOMPATIBLE");
      exit(0);
   }
   CWmatrix<T> temp(nra,ncb,0.0);
#ifdef _OPENMP
#pragma omp parallel for \
   shared(a,b,temp)\
   private(i,j,k)\
   schedule(guided)
#endif
  for (i=0; i<nra; i++){
    for (j=0; j<ncb; j++){
      for (k=0; k<nca; k++){
        temp[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  return temp;
}

template<class T>
CWmatrix<T> operator*(double b, const CWmatrix<T> a)
{
  int i,j;
  int nra = a.nrows();
  int nca = a.ncols();

  CWmatrix<T> temp(nra,nca,0.0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i=0; i<nra; i++){
    for (j=0; j<nca; j++){
       temp[i][j] = a[i][j] * b;
    }
  }
  return temp;
}

/* Matrix x a real number */
template <class T>
CWmatrix<T> operator*(const CWmatrix<T> a, double b)
{
  int i,j;
  int nra = a.nrows();
  int nca = a.ncols();

  CWmatrix<T> temp(nra,nca,0.0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i=0; i<nra; i++){
    for (j=0; j<nca; j++){
       temp[i][j] = a[i][j] * b;
    }
  }
  return temp;
}


/*Transpose in place: A<-A^T*/
template <class T>
inline void CWmatrix<T>::trans()
{
   int i,j,temp,nel;
   T **t;
   t = mm>0 ? new T*[mm] : NULL;
   nel = mm*nn;
   if (t) t[0] = nel>0 ? new T[nel] : NULL;
   for (i=1; i< mm; i++) t[i] = t[i-1] + nn;
   // Copy v to t
   for (i=0; i<mm; i++) for (j=0; j<nn; j++) t[i][j]=v[j][i];
   // resize v to (mm,nn)
   this->resize(mm,nn);
   for (i=0; i<nn; i++) {
      for (j=0; j<mm; j++) {
         v[i][j] = t[i][j];
      }
   }
   delete[] (t[0]);
   delete[] (t);
}

template <class T>
inline void CWmatrix<T>::trans(CWmatrix<T> temp)
{
   int i,j;
   int nra = mm;
   int nca = nn;
   //temp = new CWmatrix<T>();
   temp.assign(nra,nca,0.0);
   for (i=0; i<nn; i++) {
      for (j=0; j<mm; j++) {
         temp[j][i] = v[i][j];
      }
   }
   return temp;
}


template <class T>
CWmatrix<T>::~CWmatrix()
{
  if (v != NULL) {
    delete[] (v[0]);
    delete[] (v);
  }
}

// for convinience
typedef int Int; // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer

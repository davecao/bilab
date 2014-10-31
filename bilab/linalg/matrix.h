#ifndef MATRIX_H
#define MATRIX_H

#include "common.h"

using namespace std;

// pre-declare the template 
template<typename T> class CWmatrix;

template<typename T>
CWmatrix<T> operator+ (const CWmatrix<T> a, const CWmatrix<T> b);

template<typename T>
CWmatrix<T> operator- (const CWmatrix<T> a, const CWmatrix<T> b);

template<typename T>
CWmatrix<T> operator* (const CWmatrix<T> a, const CWmatrix<T> b);

template<typename T>
CWmatrix<T> operator* (const CWmatrix<T> a, double b);

template<typename T>
CWmatrix<T> operator* (double b, const CWmatrix<T> a);


template <class T>
class CWmatrix {

private:

  int nn;
  int mm;
  T **v;

public:
  CWmatrix();
  CWmatrix(int n, int m);     // Zero-based array
  CWmatrix(int n, int m, const T &a); //Initialize to constant
  CWmatrix(int n, int m, const T *a); // Initialize to array
  CWmatrix(const CWmatrix &rhs);    // Copy constructor

  ~CWmatrix();

  //Assignment operators
  CWmatrix & operator=(const CWmatrix &rhs);
  CWmatrix & operator=(const T& a);

  // The following friend template functions could not be compiled in pysetup
  // Matrix to Matrix: operators overload
//  friend CWmatrix<T> operator+<>(const CWmatrix, const CWmatrix);
//  friend CWmatrix<T> operator-<>(const CWmatrix, const CWmatrix);
//  friend CWmatrix<T> operator*<>(const CWmatrix, const CWmatrix);

    //scalar multiplication and division
//  friend CWmatrix<T> operator*<>(const CWmatrix, double);// matrix x real
//  friend CWmatrix<T> operator*<>(double, const CWmatrix);// real x matrix

  // Declare templates above
  // Matrix to Matrix: operators overload
  friend CWmatrix<T> operator+<>(const CWmatrix<T> a, const CWmatrix<T> b);
  friend CWmatrix<T> operator-<>(const CWmatrix<T> a, const CWmatrix<T> b);
  friend CWmatrix<T> operator*<>(const CWmatrix<T> a, const CWmatrix<T> b);

    //scalar multiplication and division
  friend CWmatrix<T> operator*<>(const CWmatrix<T> a, double b);// matrix x real
  friend CWmatrix<T> operator*<>(double b, const CWmatrix<T> a);// real x matrix

  //friend bool operator==(const CWmatrix);
  //friend bool operator!() const { return iszero();};

  typedef T value_type; // make T available externally
  inline T* operator[](const int i);  //subscripting: pointer to row i
  inline const T* operator[](const int i) const;
  inline T& operator()(const int i, const int j);  //subscripting: pointer to row i
  inline const T& operator()(const int i, const int j) const;
  inline int nrows() const;
  inline int ncols() const;
  inline void trans(); // matrix transpose in place
  inline void trans(const CWmatrix);// matrix transpose
  inline void show() const;
  inline void show(string info) const;

  void resize(int newn, int newm); // resize (contents not preserved)
  void assign(int newn, int newm, const T &a); // resize and assign a constant value
  T*   getData(){return *v;};
  void CleanUp(){ this->v=NULL; }; // explicit call deconstructor

};

#endif

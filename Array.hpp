#ifndef ARRAY_H
#define ARRAY_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdio.h>
#include <vector>

template <class T>
class Array {
 public:
  Array<T>(int, int);  // Constructor
  ~Array<T>() {};  // Destructor

  // Arithmetic operators: two arrays
  Array<T> operator+(const Array<T>&) const;
  Array<T> operator-(const Array<T>&) const;
  Array<T> operator*(const Array<T>&) const;
  Array<T> operator/(const Array<T>&) const;

  // Assignment operators: two arrays
  Array<T>& operator+=(const Array<T>&);
  Array<T>& operator-=(const Array<T>&);
  Array<T>& operator*=(const Array<T>&);
  Array<T>& operator/=(const Array<T>&);

  // Arithmetic operators: array and constant
  Array<T> operator+(const T&) const;
  Array<T> operator-(const T&) const;
  Array<T> operator*(const T&) const;
  Array<T> operator/(const T&) const;

  // Assignment operators: array and constant
  Array<T>& operator+=(const T&);
  Array<T>& operator-=(const T&);
  Array<T>& operator*=(const T&);
  Array<T>& operator/=(const T&);

  // Fill an array
  Array<T>& fill(const Array<T>&);
  Array<T>& fill(const T&);

  // Data access
  T& operator()(const int&, const int&);
  const T& operator()(const int&, const int&) const;

  // Data access (1D)
  T& operator()(const int&);
  const T& operator()(const int&) const;

  void print(int, int);
  void print(int, int, const char*);
  void printsci(int, const char*);

  int get_nx() const {return nx;};
  int get_ny() const {return ny;};
  int get_n() const {return n;};

 private:
  int nx, ny, n;
  std::vector<T> data;
};

#endif

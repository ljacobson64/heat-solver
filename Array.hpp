#ifndef ARRAY_H
#define ARRAY_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdio.h>
#include <vector>

class Array {
 public:
  Array(int, int);  // Constructor
  ~Array() {};  // Destructor

  // Arithmetic operators: two arrays
  Array operator+(const Array&) const;
  Array operator-(const Array&) const;
  Array operator*(const Array&) const;
  Array operator/(const Array&) const;

  // Assignment operators: two arrays
  Array& operator+=(const Array&);
  Array& operator-=(const Array&);
  Array& operator*=(const Array&);
  Array& operator/=(const Array&);

  // Arithmetic operators: array and constant
  Array operator+(const double&) const;
  Array operator-(const double&) const;
  Array operator*(const double&) const;
  Array operator/(const double&) const;

  // Assignment operators: array and constant
  Array& operator+=(const double&);
  Array& operator-=(const double&);
  Array& operator*=(const double&);
  Array& operator/=(const double&);

  // Fill an array
  Array& fill(const Array&);
  Array& fill(const double&);

  // Data access
  double& operator()(const int&, const int&);
  const double& operator()(const int&, const int&) const;

  // Data access (1D)
  double& operator()(const int&);
  const double& operator()(const int&) const;

  void print(int, int);

  int get_nx() {return nx;};
  int get_ny() {return ny;};
  int get_n() {return n;};

 private:
  int nx, ny, n;
  std::vector<double> data;
};

#endif

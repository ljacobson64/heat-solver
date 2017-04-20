#ifndef ARRAY_H
#define ARRAY_H

#include <cassert>
#include <iostream>
#include <stdio.h>

class Array {
  public:
    Array(int, int);  // Constructor
    ~Array();  // Destructor

    Array* operator+(Array*);
    Array* operator*(Array*);
    Array* operator+=(Array*);
    Array* operator*=(Array*);

    Array* operator+(double);
    Array* operator*(double);
    Array* operator+=(double);
    Array* operator*=(double);

    void fill(Array*);
    void fill(double);

    void setVal(int, int, double);
    double getVal(int, int);
    double getVal(int);
    void print(int, int);

    int get_nx() {return nx;};
    int get_ny() {return ny;};
    int get_n() {return n;};

  private:
    int nx, ny, n;
    double* data;
};

#endif

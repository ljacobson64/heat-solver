#include "Array.hpp"

// Constructor
Array::Array(int cols, int rows) {
  assert(cols > 0 && rows > 0);
  nx = cols;
  ny = rows;
  n = nx * ny;
  data = new double[n];
  fill((double)0);
}

// Destructor
Array::~Array() { delete[] data; }

// Addition operator for two arrays
Array* Array::operator+(Array* other) {
  assert(nx == other->nx && ny == other->ny);
  Array* result;
  result = new Array(nx, ny);
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      result->data[j*nx + i] = data[j*nx + i] + other->data[j*nx + i];
  return result;
}

// Multiplication operator for two arrays
Array* Array::operator*(Array* other) {
  assert(nx == other->nx && ny == other->ny);
  Array* result;
  result = new Array(nx, ny);
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      result->data[j*nx + i] = data[j*nx + i] * other->data[j*nx + i];
  return result;
}

// Compound addition assignment operator for two arrays
Array* Array::operator+=(Array* other) {
  assert(nx == other->nx && ny == other->ny);
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      data[j*nx + i] += other->data[j*nx + i];
}

// Compound multiplication assignment operator for two arrays
Array* Array::operator*=(Array* other) {
  assert(nx == other->nx && ny == other->ny);
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      data[j*nx + i] *= other->data[j*nx + i];
}

// Addition operator for an array and a constant
Array* Array::operator+(double constant) {
  Array* result;
  result = new Array(nx, ny);
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      result->data[j*nx + i] = data[j*nx + i] + constant;
  return result;
}

// Multiplication operator for an array and a constant
Array* Array::operator*(double constant) {
  Array* result;
  result = new Array(nx, ny);
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      result->data[j*nx + i] = data[j*nx + i] * constant;
  return result;
}

// Compound addition assignment operator an array and a constant
Array* Array::operator+=(double constant) {
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      data[j*nx + i] += constant;
}

// Compound multiplication assignment operator an array and a constant
Array* Array::operator*=(double constant) {
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      data[j*nx + i] *= constant;
}

// Fill an array with data copied from another
void Array::fill(Array* other) {
  assert(nx == other->nx && ny == other->ny);
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      data[j*nx + i] = other->data[j*nx + i];
}

// Fill an array with a constant value
void Array::fill(double constant) {
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      data[j*nx + i] = constant;
}

void Array::setVal(int i, int j, double val) {
  assert(i >= 0 && i < nx && j >= 0 && j < ny);
  data[j*nx + i] = val;
}

double Array::getVal(int i, int j) {
  assert(i >= 0 && i < nx && j >= 0 && j < ny);
  return data[j*nx + i];
}

double Array::getVal(int k) {
  assert((nx == 1 || ny == 1) && k < n);
  return data[k];
}

void Array::print(int characters, int decimals) {
  for (int j = 0; j < ny; j++) {
    std::string format_string = "%" + std::to_string(characters) +
                                "." + std::to_string(decimals) + "f ";
    for (int i = 0; i < nx; i++) {
      printf(format_string.c_str(), getVal(i, j));
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

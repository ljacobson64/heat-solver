#include "Array.hpp"

// Constructor
Array::Array(int cols, int rows) {
  assert(cols > 0 && rows > 0);
  nx = cols;
  ny = rows;
  n = nx * ny;
  data.assign(n, 0);
}

// Arithmetic operators for two arrays
Array Array::operator+(const Array& other) const {
  assert(nx == other.nx && ny == other.ny);
  Array result(nx, ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 result.data.begin(), std::plus<double>());
  return result;
}

Array Array::operator-(const Array& other) const {
  assert(nx == other.nx && ny == other.ny);
  Array result(nx, ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 result.data.begin(), std::minus<double>());
  return result;
}

Array Array::operator*(const Array& other) const {
  assert(nx == other.nx && ny == other.ny);
  Array result(nx, ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 result.data.begin(), std::multiplies<double>());
  return result;
}

Array Array::operator/(const Array& other) const {
  assert(nx == other.nx && ny == other.ny);
  Array result(nx, ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 result.data.begin(), std::divides<double>());
  return result;
}

// Assignment operators for two arrays
Array& Array::operator+=(const Array& other) {
  assert(nx == other.nx && ny == other.ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 data.begin(), std::plus<double>());
}

Array& Array::operator-=(const Array& other) {
  assert(nx == other.nx && ny == other.ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 data.begin(), std::minus<double>());
}

Array& Array::operator*=(const Array& other) {
  assert(nx == other.nx && ny == other.ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 data.begin(), std::multiplies<double>());
}

Array& Array::operator/=(const Array& other) {
  assert(nx == other.nx && ny == other.ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 data.begin(), std::divides<double>());
}

// Arithmetic operators for an array and a constant
Array Array::operator+(const double& constant) const {
  Array result(nx, ny);
  std::transform(data.begin(), data.end(), result.data.begin(),
                 bind2nd(std::plus<double>(), constant));
  return result;
}

Array Array::operator-(const double& constant) const {
  Array result(nx, ny);
  std::transform(data.begin(), data.end(), result.data.begin(),
                 bind2nd(std::minus<double>(), constant));
  return result;
}

Array Array::operator*(const double& constant) const {
  Array result(nx, ny);
  std::transform(data.begin(), data.end(), result.data.begin(),
                 bind2nd(std::multiplies<double>(), constant));
  return result;
}

Array Array::operator/(const double& constant) const {
  Array result(nx, ny);
  std::transform(data.begin(), data.end(), result.data.begin(),
                 bind2nd(std::divides<double>(), constant));
  return result;
}

// Assignment operators for an array and a constant
Array& Array::operator+=(const double& constant) {
  std::transform(data.begin(), data.end(), data.begin(),
                 bind2nd(std::plus<double>(), constant));
}

Array& Array::operator-=(const double& constant) {
  std::transform(data.begin(), data.end(), data.begin(),
                 bind2nd(std::minus<double>(), constant));
}

Array& Array::operator*=(const double& constant) {
  std::transform(data.begin(), data.end(), data.begin(),
                 bind2nd(std::multiplies<double>(), constant));
}

Array& Array::operator/=(const double& constant) {
  std::transform(data.begin(), data.end(), data.begin(),
                 bind2nd(std::divides<double>(), constant));
}

// Fill an array with data copied from another
Array& Array::fill(const Array& other) {
  assert(nx == other.nx && ny == other.ny);
  data.assign(other.data.begin(), other.data.end());
}

// Fill an array with a constant value
Array& Array::fill(const double& constant) {
  data.assign(n, constant);
}

// Data access
double& Array::operator()(const int& i, const int& j) {
  assert(i >= 0 && i < nx && j >= 0 && j < ny);
  return data[j * nx + i];
}

const double& Array::operator()(const int& i, const int& j) const {
  assert(i >= 0 && i < nx && j >= 0 && j < ny);
  return data[j * nx + i];
}

// Data access (1D)
double& Array::operator()(const int& k) {
  assert((nx == 1 || ny == 1) && k < n);
  return data[k];
}

const double& Array::operator()(const int& k) const {
  assert((nx == 1 || ny == 1) && k < n);
  return data[k];
}

void Array::print(int characters, int decimals) {
  std::string format_string = "%" + std::to_string(characters) +
                              "." + std::to_string(decimals) + "f ";
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      printf(format_string.c_str(), (*this)(i, j));
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

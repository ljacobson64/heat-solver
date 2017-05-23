#include "Array.hpp"

// Constructor
template <class T>
Array<T>::Array(int cols, int rows) {
  assert(cols > 0 && rows > 0);
  nx = cols;
  ny = rows;
  n = nx * ny;
  data.assign(n, 0);
}

// Arithmetic operators for two arrays
template <class T>
Array<T> Array<T>::operator+(const Array<T>& other) const {
  assert(nx == other.nx && ny == other.ny);
  Array<T> result(nx, ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 result.data.begin(), std::plus<T>());
  return result;
}

template <class T>
Array<T> Array<T>::operator-(const Array<T>& other) const {
  assert(nx == other.nx && ny == other.ny);
  Array<T> result(nx, ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 result.data.begin(), std::minus<T>());
  return result;
}

template <class T>
Array<T> Array<T>::operator*(const Array<T>& other) const {
  assert(nx == other.nx && ny == other.ny);
  Array<T> result(nx, ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 result.data.begin(), std::multiplies<T>());
  return result;
}

template <class T>
Array<T> Array<T>::operator/(const Array<T>& other) const {
  assert(nx == other.nx && ny == other.ny);
  Array<T> result(nx, ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 result.data.begin(), std::divides<T>());
  return result;
}

// Assignment operators for two arrays
template <class T>
Array<T>& Array<T>::operator+=(const Array<T>& other) {
  assert(nx == other.nx && ny == other.ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 data.begin(), std::plus<T>());
}

template <class T>
Array<T>& Array<T>::operator-=(const Array<T>& other) {
  assert(nx == other.nx && ny == other.ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 data.begin(), std::minus<T>());
}

template <class T>
Array<T>& Array<T>::operator*=(const Array<T>& other) {
  assert(nx == other.nx && ny == other.ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 data.begin(), std::multiplies<T>());
}

template <class T>
Array<T>& Array<T>::operator/=(const Array<T>& other) {
  assert(nx == other.nx && ny == other.ny);
  std::transform(data.begin(), data.end(), other.data.begin(),
                 data.begin(), std::divides<T>());
}

// Arithmetic operators for an array and a constant
template <class T>
Array<T> Array<T>::operator+(const T& constant) const {
  Array<T> result(nx, ny);
  std::transform(data.begin(), data.end(), result.data.begin(),
                 bind2nd(std::plus<T>(), constant));
  return result;
}

template <class T>
Array<T> Array<T>::operator-(const T& constant) const {
  Array<T> result(nx, ny);
  std::transform(data.begin(), data.end(), result.data.begin(),
                 bind2nd(std::minus<T>(), constant));
  return result;
}

template <class T>
Array<T> Array<T>::operator*(const T& constant) const {
  Array<T> result(nx, ny);
  std::transform(data.begin(), data.end(), result.data.begin(),
                 bind2nd(std::multiplies<T>(), constant));
  return result;
}

template <class T>
Array<T> Array<T>::operator/(const T& constant) const {
  Array<T> result(nx, ny);
  std::transform(data.begin(), data.end(), result.data.begin(),
                 bind2nd(std::divides<T>(), constant));
  return result;
}

// Assignment operators for an array and a constant
template <class T>
Array<T>& Array<T>::operator+=(const T& constant) {
  std::transform(data.begin(), data.end(), data.begin(),
                 bind2nd(std::plus<T>(), constant));
}

template <class T>
Array<T>& Array<T>::operator-=(const T& constant) {
  std::transform(data.begin(), data.end(), data.begin(),
                 bind2nd(std::minus<T>(), constant));
}

template <class T>
Array<T>& Array<T>::operator*=(const T& constant) {
  std::transform(data.begin(), data.end(), data.begin(),
                 bind2nd(std::multiplies<T>(), constant));
}

template <class T>
Array<T>& Array<T>::operator/=(const T& constant) {
  std::transform(data.begin(), data.end(), data.begin(),
                 bind2nd(std::divides<T>(), constant));
}

// Fill an array with data copied from another
template <class T>
Array<T>& Array<T>::fill(const Array<T>& other) {
  assert(nx == other.nx && ny == other.ny);
  data.assign(other.data.begin(), other.data.end());
}

// Fill an array with a constant value
template <class T>
Array<T>& Array<T>::fill(const T& constant) {
  data.assign(n, constant);
}

// Data access
template <class T>
T& Array<T>::operator()(const int& i, const int& j) {
  assert(i >= 0 && i < nx && j >= 0 && j < ny);
  return data[j * nx + i];
}

template <class T>
const T& Array<T>::operator()(const int& i, const int& j) const {
  assert(i >= 0 && i < nx && j >= 0 && j < ny);
  return data[j * nx + i];
}

// Data access (1D)
template <class T>
T& Array<T>::operator()(const int& k) {
  assert((nx == 1 || ny == 1) && k < n);
  return data[k];
}

template <class T>
const T& Array<T>::operator()(const int& k) const {
  assert((nx == 1 || ny == 1) && k < n);
  return data[k];
}

template <class T>
void Array<T>::print(int characters, int decimals) {
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

template <class T>
void Array<T>::print(int characters, int decimals, const char* fname) {
  std::string format_string = "%" + std::to_string(characters) +
                              "." + std::to_string(decimals) + "f ";
  FILE* outFile;
  outFile = fopen(fname, "w");
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      fprintf(outFile, format_string.c_str(), (*this)(i, j));
    }
    fprintf(outFile, "\n");
  }
  fprintf(outFile, "\n");
}

template <class T>
void Array<T>::printsci(int decimals) {
  std::string format_string = "%" + std::to_string(decimals + 6) +
                              "." + std::to_string(decimals) + "e ";
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      printf(format_string.c_str(), (*this)(i, j));
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template <class T>
void Array<T>::printsci(int decimals, const char* fname) {
  std::string format_string = "%" + std::to_string(decimals + 6) +
                              "." + std::to_string(decimals) + "e ";
  FILE* outFile;
  outFile = fopen(fname, "w");
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      fprintf(outFile, format_string.c_str(), (*this)(i, j));
    }
    fprintf(outFile, "\n");
  }
  fprintf(outFile, "\n");
}

template class Array<double>;

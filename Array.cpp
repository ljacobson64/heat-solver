#include "Array.hpp"

// Display an error message and exit
void throw_error(const std::string msg) {
  std::cerr << BOLD(FRED("ERROR: ")) << msg << std::endl;
  exit(EXIT_FAILURE);
}

// Display a warning message
void throw_warning(const std::string msg) {
  std::cerr << BOLD(FYEL("WARNING: ")) << msg << std::endl;
}

// Open a file
std::ifstream open_file(const std::string fname) {
  std::ifstream infile(fname);
  if (infile.fail())
    throw_error(fname + " not found");
  return infile;
}

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

// Fill an array with data read from file
template <class T>
Array<T>& Array<T>::fill_from_file(const std::string fname) {
  std::ifstream infile = open_file(fname);
  std::string buffer;
  T temp;
  int i = 0, j = 0;
  std::getline(infile, buffer);
  while (!buffer.empty()) {
    if (j == ny)
      throw_error("Read " + std::to_string(j + 1) +
                  " lines, expected " + std::to_string(ny));
    std::stringstream iss(buffer);
    i = 0;
    while (iss >> temp) {
      if (i == nx)
        throw_error("Read " + std::to_string(i + 1) +
                    " values on line, expected " + std::to_string(nx));
      data[j * nx + i] = temp;
      i++;
    }
    if (i < nx)
      throw_error("Read " + std::to_string(i) +
                  " values on line, expected " + std::to_string(nx));
    std::getline(infile, buffer);
    j++;
  }
  if (j < ny)
    throw_error("Read " + std::to_string(j) +
                " lines, expected " + std::to_string(ny));
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
                              "." + std::to_string(decimals) + "f";
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx - 1; i++)
      printf((format_string + " ").c_str(), (*this)(i, j));
    printf(format_string.c_str(), (*this)(nx - 1, j));
    std::cout << std::endl;
  }
}

template <class T>
void Array<T>::print(int characters, int decimals, const std::string fname) {
  std::string format_string = "%" + std::to_string(characters) +
                              "." + std::to_string(decimals) + "f";
  FILE* outFile;
  outFile = fopen(fname.c_str(), "w");
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx - 1; i++)
      fprintf(outFile, (format_string + " ").c_str(), (*this)(i, j));
    fprintf(outFile, format_string.c_str(), (*this)(nx - 1, j));
    fprintf(outFile, "\n");
  }
}

template <class T>
void Array<T>::printsci(int decimals) {
  std::string format_string = "%" + std::to_string(decimals + 6) +
                              "." + std::to_string(decimals) + "e";
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx - 1; i++)
      printf((format_string + " ").c_str(), (*this)(i, j));
    printf(format_string.c_str(), (*this)(nx - 1, j));
    std::cout << std::endl;
  }
}

template <class T>
void Array<T>::printsci(int decimals, const std::string fname) {
  std::string format_string = "%" + std::to_string(decimals + 6) +
                              "." + std::to_string(decimals) + "e";
  FILE* outFile;
  outFile = fopen(fname.c_str(), "w");
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx - 1; i++)
      fprintf(outFile, (format_string + " ").c_str(), (*this)(i, j));
    fprintf(outFile, format_string.c_str(), (*this)(nx - 1, j));
    fprintf(outFile, "\n");
  }
}

template class Array<double>;

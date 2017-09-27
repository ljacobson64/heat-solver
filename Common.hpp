#ifndef COMMON_H
#define COMMON_H

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <math.h>

// Output color and formatting
#define RESET "\x1B[0m"
#define FRED(x) "\x1B[31m" x RESET
#define FYEL(x) "\x1B[33m" x RESET
#define BOLD(x) "\x1B[1m" x RESET

void throw_error(const std::string);
void throw_warning(const std::string);

std::ifstream open_file(std::string);

#endif

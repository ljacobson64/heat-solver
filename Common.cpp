#include "Common.hpp"

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

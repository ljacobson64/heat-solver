#include "Array.hpp"
#include "BoundaryCond.hpp"
#include "Common.hpp"

void print_usage() {
  std::cout << "Usage: ./heat [-d <data_dir>]" << std::endl;
  exit(EXIT_SUCCESS);
}

void read_params_from_file(const std::string fname, double& Lx, double& Ly,
                           int& nx, int& ny) {
  std::ifstream infile = open_file(fname);
  std::string buffer;

  // Length in x-direction [m]
  try {
    std::getline(infile, buffer);
    Lx = std::stod(buffer);
    if (Lx <= 0)
      throw 1;
  } catch (...) {
    throw_error("Lx = " + buffer + ", must be a positive number");
  }
  // Length in y-direction [m]
  try {
    std::getline(infile, buffer);
    Ly = std::stod(buffer);
    if (Ly <= 0)
      throw 1;
  } catch (...) {
    throw_error("Ly = " + buffer + ", must be a positive number");
  }
  // Number of spatial regions in x-direction
  try {
    std::getline(infile, buffer);
    nx = std::stoi(buffer);
    if (nx <= 0 || std::stod(buffer) != (double)nx)
      throw 1;
  } catch (...) {
    throw_error("nx = " + buffer + ", must be a positive integer");
  }
  // Number of spatial regions in y-direction
  try {
    std::getline(infile, buffer);
    ny = std::stoi(buffer);
    if (ny <= 0 || std::stod(buffer) != (double)ny)
      throw 1;
  } catch (...) {
    throw_error("ny = " + buffer + ", must be a positive integer");
  }

  infile.close();
}

void solve_fourier_2D(const Array<double> x, const Array<double> y,
                      const Array<double> k, const Array<double> Q,
                      BoundaryConds BC, Array<double>* T,
                      int max_it, double tol) {
  int nx = x.get_nx() - 1;
  int ny = y.get_ny() - 1;

  // dx: on intervals (1D)
  Array<double> dx(nx, 1);
  for (int i = 0; i < nx; i++)
    dx(i) = x(i + 1, 0) - x(i, 0);

  // dy: on intervals (1D)
  Array<double> dy(1, ny);
  for (int j = 0; j < ny; j++)
    dy(j) = y(0, j + 1) - y(0, j);

  // kdy: on dx intervals
  Array<double> kdy(nx, ny + 1);
  for (int j = 1; j < ny; j++)
    for (int i = 0; i < nx; i++)  // Internal
      kdy(i, j) = 0.5 * (k(i, j - 1) * dy(j - 1) + k(i, j) * dy(j));
  for (int i = 0; i < nx; i++) {  // Bottom, top
    kdy(i, 0) = 0.5 * k(i, 0) * dy(0);
    kdy(i, ny) = 0.5 * k(i, ny - 1) * dy(ny - 1);
  }

  // kdx: on dy intervals
  Array<double> kdx(nx + 1, ny);
  for (int j = 0; j < ny; j++)
    for (int i = 1; i < nx; i++)  // Internal
      kdx(i, j) = 0.5 * (k(i - 1, j) * dx(i - 1) + k(i, j) * dx(i));
  for (int j = 0; j < ny; j++) {  // Left, right
    kdx(0, j) = 0.5 * k(0, j) * dx(0);
    kdx(nx, j) = 0.5 * k(nx - 1, j) * dx(nx - 1);
  }

  // kdy/dx: on dx intervals
  Array<double> kdydx(nx, ny + 1);
  for (int j = 0; j < ny + 1; j++)
    for (int i = 0; i < nx; i++)
      kdydx(i, j) = kdy(i, j) / dx(i);

  // kdx/dy: on dy intervals
  Array<double> kdxdy(nx + 1, ny);
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx + 1; i++)
      kdxdy(i, j) = kdx(i, j) / dy(j);

  // Denominator
  Array<double> denom(nx + 1, ny + 1);
  for (int j = 1; j < ny; j++) {
    // Internal nodes
    for (int i = 1; i < nx; i++) {
      denom(i, j) = kdydx(i - 1, j) + kdydx(i, j) +
                    kdxdy(i, j - 1) + kdxdy(i, j);
    }
    // Left
    if (BC.left.get_type() == BC_FIXED_TEMP)
      denom(0, j) = 0.;
    else if (BC.left.get_type() == BC_CONVECTIVE)
      denom(0, j) = BC.left.get_h() * dy(j) + kdydx(0, j) +
                    kdxdy(0, j - 1) + kdxdy(0, j);
    // Right
    if (BC.right.get_type() == BC_FIXED_TEMP)
      denom(nx, j) = 0.;
    else if (BC.right.get_type() == BC_CONVECTIVE)
      denom(nx, j) = kdydx(nx - 1, j) + BC.right.get_h() * dy(j) +
                     kdxdy(nx, j - 1) + kdxdy(nx, j);
  }
  for (int i = 1; i < nx; i++) {
    // Bottom
    if (BC.bottom.get_type() == BC_FIXED_TEMP)
      denom(i, 0) = 0.;
    else if (BC.bottom.get_type() == BC_CONVECTIVE)
      denom(i, 0) = kdydx(i - 1, 0) + kdydx(i, 0) +
                    BC.bottom.get_h() * dx(i) + kdxdy(i, 0);
    // Top
    if (BC.top.get_type() == BC_FIXED_TEMP)
      denom(i, ny) = 0.;
    else if (BC.top.get_type() == BC_CONVECTIVE)
      denom(i, ny) = kdydx(i - 1, ny) + kdydx(i, ny) +
                     kdxdy(i, ny - 1) + BC.top.get_h() * dx(i);
  }
  // Corners
  if (BC.left.get_type() == BC_CONVECTIVE && BC.bottom.get_type() == BC_CONVECTIVE)
    denom(0, 0) = 0.5 * BC.left.get_h() * dy(0) + kdydx(0, 0) +
                  0.5 * BC.bottom.get_h() * dx(0) + kdxdy(0, 0);
  if (BC.right.get_type() == BC_CONVECTIVE && BC.bottom.get_type() == BC_CONVECTIVE)
    denom(nx, 0) = kdydx(nx - 1, 0) + 0.5 * BC.right.get_h() * dy(0) +
                   0.5 * BC.bottom.get_h() * dx(nx - 1) + kdxdy(nx, 0);
  if (BC.left.get_type() == BC_CONVECTIVE && BC.top.get_type() == BC_CONVECTIVE)
    denom(0, ny) = 0.5 * BC.left.get_h() * dy(ny - 1) + kdydx(0, ny) +
                   kdxdy(0, ny - 1) + 0.5 * BC.top.get_h() * dx(0);
  if (BC.right.get_type() == BC_CONVECTIVE && BC.top.get_type() == BC_CONVECTIVE)
    denom(nx, ny) = kdydx(nx - 1, ny) + 0.5 * BC.right.get_h() * dy(ny - 1) +
                    kdxdy(nx, ny - 1) + 0.5 * BC.top.get_h() * dx(nx - 1);

  // Heat source
  Array<double> Q_gen(nx + 1, ny + 1);
  for (int j = 1; j < ny; j++) {
    // Internal nodes
    for (int i = 1; i < nx; i++)
      Q_gen(i, j) = 0.25 * Q(i, j) * (dx(i - 1) + dx(i)) * (dy(j - 1) + dy(j));
    // Left, right
    Q_gen(0, j) = 0.25 * Q(0, j) * dx(0) * (dy(j - 1) + dy(j));
    Q_gen(nx, j) = 0.25 * Q(nx, j) * dx(nx - 1) * (dy(j - 1) + dy(j));
  }
  for (int i = 1; i < nx; i++) {
    // Bottom, top
    Q_gen(i, 0) = 0.25 * Q(i, 0) * (dx(i - 1) + dx(i)) * dy(0);
    Q_gen(i, ny) = 0.25 * Q(i, ny) * (dx(i - 1) + dx(i)) * dy(ny - 1);
  }
  // Corners
  Q_gen(0, 0) = 0.25 * Q(0, 0) * dx(0) * dy(0);
  Q_gen(nx, 0) = 0.25 * Q(nx - 1, 0) * dx(nx - 1) * dy(0);
  Q_gen(0, ny) = 0.25 * Q(0, ny - 1) * dx(0) * dy(ny - 1);
  Q_gen(nx, ny) = 0.25 * Q(nx - 1, ny - 1) * dx(nx - 1) * dy(ny - 1);

  Array<double> T_old(nx + 1, ny + 1);

  int it = 0;
  double max_dif = 1e100;
  double omega = 1.6;

  // Set fixed temperatures as appropriate
  if (BC.left.get_type() == BC_FIXED_TEMP)
    for (int j = 0; j < ny + 1; j++)
      (*T)(0, j) = BC.left.get_T_fixed();
  if (BC.right.get_type() == BC_FIXED_TEMP)
    for (int j = 0; j < ny + 1; j++)
      (*T)(nx, j) = BC.right.get_T_fixed();
  if (BC.bottom.get_type() == BC_FIXED_TEMP)
    for (int i = 0; i < nx + 1; i++)
      (*T)(i, 0) = BC.bottom.get_T_fixed();
  if (BC.top.get_type() == BC_FIXED_TEMP)
    for (int i = 0; i < nx + 1; i++)
      (*T)(i, ny) = BC.top.get_T_fixed();
  // If adjacent sides have fixed temperatures, take the mean of them at the corner
  if (BC.left.get_type() == BC_FIXED_TEMP && BC.bottom.get_type() == BC_FIXED_TEMP)
    (*T)(0, 0) = 0.5 * (BC.left.get_T_fixed() + BC.bottom.get_T_fixed());
  if (BC.right.get_type() == BC_FIXED_TEMP && BC.bottom.get_type() == BC_FIXED_TEMP)
    (*T)(nx, 0) = 0.5 * (BC.right.get_T_fixed() + BC.bottom.get_T_fixed());
  if (BC.left.get_type() == BC_FIXED_TEMP && BC.top.get_type() == BC_FIXED_TEMP)
    (*T)(0, ny) = 0.5 * (BC.left.get_T_fixed() + BC.top.get_T_fixed());
  if (BC.right.get_type() == BC_FIXED_TEMP && BC.top.get_type() == BC_FIXED_TEMP) {
    (*T)(nx, ny) = 0.5 * (BC.right.get_T_fixed() + BC.top.get_T_fixed());
  }

  while (max_dif > tol && it < max_it) {
    it++;

    // Old temperature values
    T_old.fill(*T);

    // Bottom-left corner
    if (BC.bottom.get_type() == BC_CONVECTIVE && BC.left.get_type() == BC_CONVECTIVE) {
      (*T)(0, 0) = (0.5 * BC.left.get_h() * BC.left.get_T_inf() * dy(0) +
                    kdydx(0, 0) * (*T)(1, 0) +
                    0.5 * BC.bottom.get_h() * BC.bottom.get_T_inf() * dx(0) +
                    kdxdy(0, 0) * (*T)(0, 1) +
                    Q_gen(0, 0)) / denom(0, 0);
    }
    // Bottom boundary
    if (BC.bottom.get_type() == BC_CONVECTIVE) {
      for (int i = 1; i < nx; i++) {
        (*T)(i, 0) = (kdydx(i - 1, 0) * (*T)(i - 1, 0) +
                      kdydx(i, 0) * (*T)(i + 1, 0) +
                      BC.bottom.get_h() * BC.bottom.get_T_inf() * dx(i) +
                      kdxdy(i, 0) * (*T)(i, 1) +
                      Q_gen(i, 0)) / denom(i, 0);
      }
    }
    // Bottom-right corner
    if (BC.bottom.get_type() == BC_CONVECTIVE && BC.right.get_type() == BC_CONVECTIVE) {
      (*T)(nx, 0) = (kdydx(nx - 1, 0) * (*T)(nx - 1, 0) +
                     0.5 * BC.right.get_h() * BC.right.get_T_inf() * dy(0) +
                     0.5 * BC.bottom.get_h() * BC.bottom.get_T_inf() * dx(nx - 1) +
                     kdxdy(nx, 0) * (*T)(nx, 1) +
                     Q_gen(nx, 0)) / denom(nx, 0);
    }
    for (int j = 1; j < ny; j++) {
      // Left boundary
      if (BC.left.get_type() == BC_CONVECTIVE) {
        (*T)(0, j) = (BC.left.get_h() * BC.left.get_T_inf() * dy(j) +
                      kdydx(0, j) * (*T)(1, j) +
                      kdxdy(0, j - 1) * (*T)(0, j - 1) +
                      kdxdy(0, j) * (*T)(0, j + 1) +
                      Q_gen(0, j)) / denom(0, j);
      }
      // Internal nodes
      for (int i = 1; i < nx; i++) {
        (*T)(i, j) = (kdydx(i - 1, j) * (*T)(i - 1, j) +
                      kdydx(i, j) * (*T)(i + 1, j) +
                      kdxdy(i, j - 1) * (*T)(i, j - 1) +
                      kdxdy(i, j) * (*T)(i, j + 1) +
                      Q_gen(i, j)) / denom(i, j);
      }
      // Right boundary
      if (BC.right.get_type() == BC_CONVECTIVE) {
        (*T)(nx, j) = (kdydx(nx - 1, j) * (*T)(nx - 1, j) +
                       BC.right.get_h() * BC.right.get_T_inf() * dy(j) +
                       kdxdy(nx, j - 1) * (*T)(nx, j - 1) +
                       kdxdy(nx, j) * (*T)(nx, j + 1) +
                       Q_gen(nx, j)) / denom(nx, j);
      }
    }
    // Top-left corner
    if (BC.top.get_type() == BC_CONVECTIVE && BC.left.get_type() == BC_CONVECTIVE) {
      (*T)(0, ny) = (0.5 * BC.left.get_h() * BC.left.get_T_inf() * dy(ny - 1) +
                     kdydx(0, ny) * (*T)(0 + 1, ny) +
                     kdxdy(0, ny - 1) * (*T)(0, ny - 1) +
                     0.5 * BC.top.get_h() * BC.top.get_T_inf() * dx(0) +
                     Q_gen(0, ny)) / denom(0, ny);
    }
    // Top boundary
    if (BC.top.get_type() == BC_CONVECTIVE) {
      for (int i = 1; i < nx; i++) {
        (*T)(i, ny) = (kdydx(i - 1, ny) * (*T)(i - 1, ny) +
                       kdydx(i, ny) * (*T)(i + 1, ny) +
                       kdxdy(i, ny - 1) * (*T)(i, ny - 1) +
                       BC.top.get_h() * BC.top.get_T_inf() * dx(i) +
                       Q_gen(i, ny)) / denom(i, ny);
      }
    }
    // Top-right corner
    if (BC.top.get_type() == BC_CONVECTIVE && BC.right.get_type() == BC_CONVECTIVE) {
      (*T)(nx, ny) = (kdydx(nx - 1, ny) * (*T)(nx - 1, ny) +
                      0.5 * BC.right.get_h() * BC.right.get_T_inf() * dy(ny - 1) +
                      kdxdy(nx, ny - 1) * (*T)(nx, ny - 1) +
                      0.5 * BC.top.get_h() * BC.top.get_T_inf() * dx(nx - 1) +
                      Q_gen(nx, ny)) / denom(nx, ny);
    }

    // Successive over-relaxation
    *T *= omega;
    *T += T_old * (1 - omega);

    // Convergence criterion
    double dif = 0;
    max_dif = 0;
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        dif = fabs((*T)(i, j) - T_old(i, j)) / T_old(i, j);
        if (dif > max_dif)
          max_dif = dif;
      }
    }
  }

  std::cout << "Number of iterations: " << it << ", difference: " << max_dif << std::endl;
}

int main(int argc, char** argv) {
  // Parse arguments
  std::string data_dir;
  if (argc == 1) {
    data_dir = "data";
  } else if (argc == 3) {
    if (std::strcmp(argv[1], "-d") == 0)
      data_dir = argv[2];
    else
      print_usage();
  } else {
    print_usage();
  }

  // Parameters
  double Lx, Ly;
  int nx, ny;
  read_params_from_file(data_dir + "/params.txt", Lx, Ly, nx, ny);

  int n = nx * ny;      // Total number of spatial regions
  double dx = Lx / nx;  // x-direction interval
  double dy = Ly / ny;  // y-direction interval
  double A = dx * dy;   // Area of spatial region

  // x-direction nodes [m]
  Array<double> x(nx + 1, 1);
  for (int i = 0; i < nx + 1; i++)
    x(i) = i * dx;

  // y-direction nodes [m]
  Array<double> y(1, ny + 1);
  for (int j = 0; j < ny + 1; j++)
    y(j) = j * dy;

  // Thermal conductivity [W/m-K]
  Array<double> k(nx, ny);
  k.fill_from_file(data_dir + "/k.txt");

  // Volumetric source (forward: [W/m^3], adjoint: [-])
  Array<double> Q(nx + 1, ny + 1);
  Q.fill_from_file(data_dir + "/Q.txt");

  // Boundary conditions
  BoundaryConds BC;
  BC.left.set_type(BC_CONVECTIVE, 10., 0.);
  BC.right.set_type(BC_CONVECTIVE, 10., 0.);
  BC.bottom.set_type(BC_CONVECTIVE, 10., 0.);
  BC.top.set_type(BC_CONVECTIVE, 10., 0.);
  //BC.left.set_type(BC_FIXED_TEMP, 2000.);
  //BC.right.set_type(BC_FIXED_TEMP, 2000.);
  //BC.bottom.set_type(BC_FIXED_TEMP, 2000.);
  //BC.top.set_type(BC_FIXED_TEMP, 2000.);

  // Solve
  Array<double> T(nx + 1, ny + 1);
  std::chrono::high_resolution_clock::time_point t1, t2;
  t1 = std::chrono::high_resolution_clock::now();
  solve_fourier_2D(x, y, k, Q, BC, &T, 100000, 1.e-8);
  t2 = std::chrono::high_resolution_clock::now();
  int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  std::cout << "Elapsed time: " << elapsed << " ms" << std::endl;

  // Write some arrays to file
  //k.printsci(4, data_dir + "/k.txt");
  //Q.printsci(4, data_dir + "/Q.txt");
  T.printsci(4, data_dir + "/T.txt");

  return 0;
}

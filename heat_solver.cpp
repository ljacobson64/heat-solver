#include "Array.hpp"
#include "BoundaryCond.hpp"
#include "Common.hpp"

void print_usage() {
  std::cout << "Usage: ./heat [-d <data_dir>]" << std::endl;
  exit(EXIT_SUCCESS);
}

void read_params_from_file(const std::string fname, double& Lx, double& Ly,
                           int& nx, int& ny, BoundaryConds& BCs) {
  std::ifstream infile = open_file(fname);
  std::string buffer, buffer2;

  // Length in x-direction [m]
  try {
    std::getline(infile, buffer);
    Lx = std::stod(buffer);
    if (Lx <= 0)
      throw 1;
  } catch (...) {
    throw_error("Invalid value for Lx. Expected positive number, received " + buffer);
  }
  // Length in y-direction [m]
  try {
    std::getline(infile, buffer);
    Ly = std::stod(buffer);
    if (Ly <= 0)
      throw 1;
  } catch (...) {
    throw_error("Invalid value for Ly. Expected positive number, received " + buffer);
  }
  // Number of spatial regions in x-direction
  try {
    std::getline(infile, buffer);
    nx = std::stoi(buffer);
    if (nx <= 0 || std::stod(buffer) != (double)nx)
      throw 1;
  } catch (...) {
    throw_error("Invalid value for nx. Expected positive integer, received " + buffer);
  }
  // Number of spatial regions in y-direction
  try {
    std::getline(infile, buffer);
    ny = std::stoi(buffer);
    if (ny <= 0 || std::stod(buffer) != (double)ny)
      throw 1;
  } catch (...) {
    throw_error("Invalid value for ny. Expected positive integer, received " + buffer);
  }

  // Expect blank line
  std::getline(infile, buffer);
  if (!buffer.empty())
    throw_error("Expected blank line");

  // Boundary conditions
  BCs.resize(BC_SIDE_MAX - BC_SIDE_MIN + 1);
  for (int i = BC_SIDE_MIN; i <= BC_SIDE_MAX; i++) {
    std::getline(infile, buffer);
    std::istringstream iss(buffer);
    std::vector<std::string> tokens;
    copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter(tokens));

    int num_tokens = tokens.size();
    int type = -1;
    if (num_tokens == 0) {
      throw_error("No BC tokens found for side " + std::to_string(i));
    } else {
      try {
        type = std::stoi(tokens[0]);
        if (type < BC_TYPE_MIN || type > BC_TYPE_MAX || std::stod(tokens[0]) != (double)type)
          throw 1;
      } catch (...) {
        throw_error("Invalid BC type for side " + std::to_string(i) +
                    ". Expected integer between 0 and 3, received " + tokens[0]);
      }
    }

    if (type == BC_ADIABATIC && num_tokens == 1) {
      BCs[i].set_type(type);
    } else if ((type == BC_FIXED_TEMP || type == BC_FIXED_FLUX) && num_tokens == 2) {
      double p1 = std::stod(tokens[1]);
      BCs[i].set_type(type, p1);
    } else if (type == BC_CONVECTIVE && num_tokens == 3) {
      double p1 = std::stod(tokens[1]);
      double p2 = std::stod(tokens[2]);
      BCs[i].set_type(type, p1, p2);
    } else {
      throw_error("Invalid number of BC tokens for side " + std::to_string(i));
    }
  }

  infile.close();
}

void solve_fourier_2D(const Array<double>& x, const Array<double>& y,
                      const Array<double>& k, const Array<double>& Q,
                      const BoundaryConds& BCs, Array<double>& T,
                      const int max_it, const double tol) {
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
    if (BCs[BC_LEFT].get_type() == BC_FIXED_TEMP)
      denom(0, j) = 0.;
    else if (BCs[BC_LEFT].get_type() == BC_CONVECTIVE)
      denom(0, j) = BCs[BC_LEFT].get_h() * dy(j) + kdydx(0, j) +
                    kdxdy(0, j - 1) + kdxdy(0, j);
    // Right
    if (BCs[BC_RIGHT].get_type() == BC_FIXED_TEMP)
      denom(nx, j) = 0.;
    else if (BCs[BC_RIGHT].get_type() == BC_CONVECTIVE)
      denom(nx, j) = kdydx(nx - 1, j) + BCs[BC_RIGHT].get_h() * dy(j) +
                     kdxdy(nx, j - 1) + kdxdy(nx, j);
  }
  for (int i = 1; i < nx; i++) {
    // Bottom
    if (BCs[BC_BOTTOM].get_type() == BC_FIXED_TEMP)
      denom(i, 0) = 0.;
    else if (BCs[BC_BOTTOM].get_type() == BC_CONVECTIVE)
      denom(i, 0) = kdydx(i - 1, 0) + kdydx(i, 0) +
                    BCs[BC_BOTTOM].get_h() * dx(i) + kdxdy(i, 0);
    // Top
    if (BCs[BC_TOP].get_type() == BC_FIXED_TEMP)
      denom(i, ny) = 0.;
    else if (BCs[BC_TOP].get_type() == BC_CONVECTIVE)
      denom(i, ny) = kdydx(i - 1, ny) + kdydx(i, ny) +
                     kdxdy(i, ny - 1) + BCs[BC_TOP].get_h() * dx(i);
  }
  // Corners
  if (BCs[BC_LEFT].get_type() == BC_CONVECTIVE && BCs[BC_BOTTOM].get_type() == BC_CONVECTIVE)
    denom(0, 0) = 0.5 * BCs[BC_LEFT].get_h() * dy(0) + kdydx(0, 0) +
                  0.5 * BCs[BC_BOTTOM].get_h() * dx(0) + kdxdy(0, 0);
  if (BCs[BC_RIGHT].get_type() == BC_CONVECTIVE && BCs[BC_BOTTOM].get_type() == BC_CONVECTIVE)
    denom(nx, 0) = kdydx(nx - 1, 0) + 0.5 * BCs[BC_RIGHT].get_h() * dy(0) +
                   0.5 * BCs[BC_BOTTOM].get_h() * dx(nx - 1) + kdxdy(nx, 0);
  if (BCs[BC_LEFT].get_type() == BC_CONVECTIVE && BCs[BC_TOP].get_type() == BC_CONVECTIVE)
    denom(0, ny) = 0.5 * BCs[BC_LEFT].get_h() * dy(ny - 1) + kdydx(0, ny) +
                   kdxdy(0, ny - 1) + 0.5 * BCs[BC_TOP].get_h() * dx(0);
  if (BCs[BC_RIGHT].get_type() == BC_CONVECTIVE && BCs[BC_TOP].get_type() == BC_CONVECTIVE)
    denom(nx, ny) = kdydx(nx - 1, ny) + 0.5 * BCs[BC_RIGHT].get_h() * dy(ny - 1) +
                    kdxdy(nx, ny - 1) + 0.5 * BCs[BC_TOP].get_h() * dx(nx - 1);

  // Heat source
  Array<double> Q_gen(nx + 1, ny + 1);
  for (int j = 1; j < ny; j++) {
    // Internal nodes
    for (int i = 1; i < nx; i++)
      Q_gen(i, j) = 0.25 * Q(i, j) * (dx(i - 1) + dx(i)) * (dy(j - 1) + dy(j));
    // Left, right
    Q_gen(0, j) = 0.5 * Q(0, j) * dx(0) * (dy(j - 1) + dy(j));
    Q_gen(nx, j) = 0.5 * Q(nx, j) * dx(nx - 1) * (dy(j - 1) + dy(j));
  }
  for (int i = 1; i < nx; i++) {
    // Bottom, top
    Q_gen(i, 0) = 0.5 * Q(i, 0) * (dx(i - 1) + dx(i)) * dy(0);
    Q_gen(i, ny) = 0.5 * Q(i, ny) * (dx(i - 1) + dx(i)) * dy(ny - 1);
  }
  // Corners
  Q_gen(0, 0) = Q(0, 0) * dx(0) * dy(0);
  Q_gen(nx, 0) = Q(nx, 0) * dx(nx - 1) * dy(0);
  Q_gen(0, ny) = Q(0, ny) * dx(0) * dy(ny - 1);
  Q_gen(nx, ny) = Q(nx, ny) * dx(nx - 1) * dy(ny - 1);

  Array<double> T_old(nx + 1, ny + 1);

  int it = 0;
  double max_dif = 1e100;
  double omega = 1.6;

  // Set fixed temperatures as appropriate
  if (BCs[BC_LEFT].get_type() == BC_FIXED_TEMP)
    for (int j = 0; j < ny + 1; j++)
      T(0, j) = BCs[BC_LEFT].get_T_fixed();
  if (BCs[BC_RIGHT].get_type() == BC_FIXED_TEMP)
    for (int j = 0; j < ny + 1; j++)
      T(nx, j) = BCs[BC_RIGHT].get_T_fixed();
  if (BCs[BC_BOTTOM].get_type() == BC_FIXED_TEMP)
    for (int i = 0; i < nx + 1; i++)
      T(i, 0) = BCs[BC_BOTTOM].get_T_fixed();
  if (BCs[BC_TOP].get_type() == BC_FIXED_TEMP)
    for (int i = 0; i < nx + 1; i++)
      T(i, ny) = BCs[BC_TOP].get_T_fixed();
  // If adjacent sides have fixed temperatures, take the mean of them at the corner
  if (BCs[BC_LEFT].get_type() == BC_FIXED_TEMP && BCs[BC_BOTTOM].get_type() == BC_FIXED_TEMP)
    T(0, 0) = 0.5 * (BCs[BC_LEFT].get_T_fixed() + BCs[BC_BOTTOM].get_T_fixed());
  if (BCs[BC_RIGHT].get_type() == BC_FIXED_TEMP && BCs[BC_BOTTOM].get_type() == BC_FIXED_TEMP)
    T(nx, 0) = 0.5 * (BCs[BC_RIGHT].get_T_fixed() + BCs[BC_BOTTOM].get_T_fixed());
  if (BCs[BC_LEFT].get_type() == BC_FIXED_TEMP && BCs[BC_TOP].get_type() == BC_FIXED_TEMP)
    T(0, ny) = 0.5 * (BCs[BC_LEFT].get_T_fixed() + BCs[BC_TOP].get_T_fixed());
  if (BCs[BC_RIGHT].get_type() == BC_FIXED_TEMP && BCs[BC_TOP].get_type() == BC_FIXED_TEMP)
    T(nx, ny) = 0.5 * (BCs[BC_RIGHT].get_T_fixed() + BCs[BC_TOP].get_T_fixed());

  while (max_dif > tol && it < max_it) {
    it++;

    // Old temperature values
    T_old.fill(T);

    // Bottom-left corner
    if (BCs[BC_BOTTOM].get_type() == BC_CONVECTIVE && BCs[BC_LEFT].get_type() == BC_CONVECTIVE) {
      T(0, 0) = (0.5 * BCs[BC_LEFT].get_h() * BCs[BC_LEFT].get_T_inf() * dy(0) +
                 kdydx(0, 0) * T(1, 0) +
                 0.5 * BCs[BC_BOTTOM].get_h() * BCs[BC_BOTTOM].get_T_inf() * dx(0) +
                 kdxdy(0, 0) * T(0, 1) +
                 Q_gen(0, 0)) / denom(0, 0);
    }
    // Bottom boundary
    if (BCs[BC_BOTTOM].get_type() == BC_CONVECTIVE) {
      for (int i = 1; i < nx; i++) {
        T(i, 0) = (kdydx(i - 1, 0) * T(i - 1, 0) +
                   kdydx(i, 0) * T(i + 1, 0) +
                   BCs[BC_BOTTOM].get_h() * BCs[BC_BOTTOM].get_T_inf() * dx(i) +
                   kdxdy(i, 0) * T(i, 1) +
                   Q_gen(i, 0)) / denom(i, 0);
      }
    }
    // Bottom-right corner
    if (BCs[BC_BOTTOM].get_type() == BC_CONVECTIVE && BCs[BC_RIGHT].get_type() == BC_CONVECTIVE) {
      T(nx, 0) = (kdydx(nx - 1, 0) * T(nx - 1, 0) +
                  0.5 * BCs[BC_RIGHT].get_h() * BCs[BC_RIGHT].get_T_inf() * dy(0) +
                  0.5 * BCs[BC_BOTTOM].get_h() * BCs[BC_BOTTOM].get_T_inf() * dx(nx - 1) +
                  kdxdy(nx, 0) * T(nx, 1) +
                  Q_gen(nx, 0)) / denom(nx, 0);
    }
    for (int j = 1; j < ny; j++) {
      // Left boundary
      if (BCs[BC_LEFT].get_type() == BC_CONVECTIVE) {
        T(0, j) = (BCs[BC_LEFT].get_h() * BCs[BC_LEFT].get_T_inf() * dy(j) +
                   kdydx(0, j) * T(1, j) +
                   kdxdy(0, j - 1) * T(0, j - 1) +
                   kdxdy(0, j) * T(0, j + 1) +
                   Q_gen(0, j)) / denom(0, j);
      }
      // Internal nodes
      for (int i = 1; i < nx; i++) {
        T(i, j) = (kdydx(i - 1, j) * T(i - 1, j) +
                   kdydx(i, j) * T(i + 1, j) +
                   kdxdy(i, j - 1) * T(i, j - 1) +
                   kdxdy(i, j) * T(i, j + 1) +
                   Q_gen(i, j)) / denom(i, j);
      }
      // Right boundary
      if (BCs[BC_RIGHT].get_type() == BC_CONVECTIVE) {
        T(nx, j) = (kdydx(nx - 1, j) * T(nx - 1, j) +
                    BCs[BC_RIGHT].get_h() * BCs[BC_RIGHT].get_T_inf() * dy(j) +
                    kdxdy(nx, j - 1) * T(nx, j - 1) +
                    kdxdy(nx, j) * T(nx, j + 1) +
                    Q_gen(nx, j)) / denom(nx, j);
      }
    }
    // Top-left corner
    if (BCs[BC_TOP].get_type() == BC_CONVECTIVE && BCs[BC_LEFT].get_type() == BC_CONVECTIVE) {
      T(0, ny) = (0.5 * BCs[BC_LEFT].get_h() * BCs[BC_LEFT].get_T_inf() * dy(ny - 1) +
                  kdydx(0, ny) * T(0 + 1, ny) +
                  kdxdy(0, ny - 1) * T(0, ny - 1) +
                  0.5 * BCs[BC_TOP].get_h() * BCs[BC_TOP].get_T_inf() * dx(0) +
                  Q_gen(0, ny)) / denom(0, ny);
    }
    // Top boundary
    if (BCs[BC_TOP].get_type() == BC_CONVECTIVE) {
      for (int i = 1; i < nx; i++) {
        T(i, ny) = (kdydx(i - 1, ny) * T(i - 1, ny) +
                    kdydx(i, ny) * T(i + 1, ny) +
                    kdxdy(i, ny - 1) * T(i, ny - 1) +
                    BCs[BC_TOP].get_h() * BCs[BC_TOP].get_T_inf() * dx(i) +
                    Q_gen(i, ny)) / denom(i, ny);
      }
    }
    // Top-right corner
    if (BCs[BC_TOP].get_type() == BC_CONVECTIVE && BCs[BC_RIGHT].get_type() == BC_CONVECTIVE) {
      T(nx, ny) = (kdydx(nx - 1, ny) * T(nx - 1, ny) +
                   0.5 * BCs[BC_RIGHT].get_h() * BCs[BC_RIGHT].get_T_inf() * dy(ny - 1) +
                   kdxdy(nx, ny - 1) * T(nx, ny - 1) +
                   0.5 * BCs[BC_TOP].get_h() * BCs[BC_TOP].get_T_inf() * dx(nx - 1) +
                   Q_gen(nx, ny)) / denom(nx, ny);
    }

    // Successive over-relaxation
    T *= omega;
    T += T_old * (1 - omega);

    // Convergence criterion
    double dif = 0;
    max_dif = 0;
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        dif = fabs(T(i, j) - T_old(i, j)) / T_old(i, j);
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
  BoundaryConds BCs;
  read_params_from_file(data_dir + "/params.txt", Lx, Ly, nx, ny, BCs);

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

  // Solve
  Array<double> T(nx + 1, ny + 1);
  std::chrono::high_resolution_clock::time_point t1, t2;
  t1 = std::chrono::high_resolution_clock::now();
  solve_fourier_2D(x, y, k, Q, BCs, T, 100000, 1.e-8);
  t2 = std::chrono::high_resolution_clock::now();
  int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  std::cout << "Elapsed time: " << elapsed << " ms" << std::endl;

  // Write some arrays to file
  //k.printsci(4, data_dir + "/k.txt");
  //Q.printsci(4, data_dir + "/Q.txt");
  T.printsci(4, data_dir + "/T.txt");

  return 0;
}

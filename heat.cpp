#include "Array.hpp"

#include <chrono>
#include <math.h>

void solve_fourier_2D(const Array<double> x, const Array<double> y,
                      const Array<double> k, const Array<double> Q,
                      double h, double T_inf, Array<double>* T) {

  int nx = x.get_nx() - 1;
  int ny = y.get_ny() - 1;

  // Make sure input dimensions match
  assert(k.get_nx() == nx && k.get_ny() == ny);
  assert(Q.get_nx() == nx && Q.get_ny() == ny);
  assert(T->get_nx() == nx + 1 && T->get_ny() == ny + 1);

  // dx: x-direction mesh; 1D; on intervals
  Array<double> dx(nx, 1);
  for (int i = 0; i < nx; i++)
    dx(i) = x(i + 1, 0) - x(i, 0);

  // dy: y-direction mesh; 1D; on intervals
  Array<double> dy(1, ny);
  for (int j = 0; j < ny; j++)
    dy(j) = y(0, j + 1) - y(0, j);

  // kdy: k times y-direction mesh: 2D; on dx intervals
  Array<double> kdy(nx, ny + 1);
  for (int j = 1; j < ny; j++)
    for (int i = 0; i < nx; i++)  // Internal
      kdy(i, j) = 0.5 * (k(i, j - 1) * dy(j - 1) + k(i, j) * dy(j));
  for (int i = 0; i < nx; i++) {  // Bottom, top
    kdy(i, 0) = 0.5 * k(i, 0) * dy(0);
    kdy(i, ny) = 0.5 * k(i, ny - 1) * dy(ny - 1);
  }

  // kdx: k times x-direction mesh: 2D; on dy intervals
  Array<double> kdx(nx + 1, ny);
  for (int j = 0; j < ny; j++)
    for (int i = 1; i < nx; i++)  // Internal
      kdx(i, j) = 0.5 * (k(i - 1, j) * dx(i - 1) + k(i, j) * dx(i));
  for (int j = 0; j < ny; j++) {  // Left, right
    kdx(0, j) = 0.5 * k(0, j) * dx(0);
    kdx(nx, j) = 0.5 * k(nx - 1, j) * dx(nx - 1);
  }

  // kdy/dx
  Array<double> kdydx(nx, ny + 1);
  for (int j = 0; j < ny + 1; j++)
    for (int i = 0; i < nx; i++)
      kdydx(i, j) = kdy(i, j) / dx(i);

  // kdx/dy
  Array<double> kdxdy(nx + 1, ny);
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx + 1; i++)
      kdxdy(i, j) = kdx(i, j) / dy(j);

  // Denominator
  Array<double> denom(nx + 1, ny + 1);
  for (int j = 1; j < ny; j++) {
    for (int i = 1; i < nx; i++)
      denom(i, j) = kdydx(i - 1, j) + kdydx(i, j) +
                    kdxdy(i, j - 1) + kdxdy(i, j);
    denom(0, j) = h * dy(j) + kdydx(0, j) +
                  kdxdy(0, j - 1) + kdxdy(0, j);
    denom(nx, j) = kdydx(nx - 1, j) + h * dy(j) +
                   kdxdy(nx, j - 1) + kdxdy(nx, j);
  }
  for (int i = 1; i < nx; i++) {
    denom(i, 0) = kdydx(i - 1, 0) + kdydx(i, 0) +
                  h * dx(i) + kdxdy(i, 0);
    denom(i, ny) = kdydx(i - 1, ny) + kdydx(i, ny) +
                   kdxdy(i, ny - 1) + h * dx(i);
  }
  denom(0, 0) = 0.5 * (h * dy(0) + h * dx(0)) +
                kdydx(0, 0) + kdxdy(0, 0);
  denom(nx, 0) = 0.5 * (h * dy(0) + h * dx(nx - 1)) +
                 kdydx(nx - 1, 0) + kdxdy(nx, 0);
  denom(0, ny) = 0.5 * (h * dy(ny - 1) + h * dx(0)) +
                 kdydx(0, ny) + kdxdy(0, ny - 1);
  denom(nx, ny) = 0.5 * (h * dy(ny - 1) + h * dx(nx - 1)) +
                  kdydx(nx - 1, ny) + kdxdy(nx, ny - 1);

  // Heat source
  Array<double> Q_gen(nx + 1, ny + 1);
  for (int j = 1; j < ny; j++) {
    for (int i = 1; i < nx; i++) {
      Q_gen(i, j) = 0.25 * (Q(i - 1, j - 1) * dx(i - 1) * dy(j - 1) +
                            Q(i - 1, j) * dx(i - 1) * dy(j) +
                            Q(i, j - 1) * dx(i) * dy(j - 1) +
                            Q(i, j) * dx(i) * dy(j));
    }
    Q_gen(0, j) = 0.25 * (Q(0, j - 1) * dx(0) * dy(j - 1) +
                          Q(0, j) * dx(0) * dy(j));
    Q_gen(nx, j) = 0.25 * (Q(nx - 1, j - 1) * dx(nx - 1) * dy(j - 1) +
                           Q(nx - 1, j) * dx(nx - 1) * dy(j));
  }
  for (int i = 1; i < nx; i++) {
    Q_gen(i, 0) = 0.25 * (Q(i - 1, 0) * dx(i - 1) * dy(0) +
                          Q(i, 0) * dx(i) * dy(0));
    Q_gen(i, ny) = 0.25 * (Q(i - 1, ny - 1) * dx(i - 1) * dy(ny - 1) +
                           Q(i, ny - 1) * dx(i) * dy(ny - 1));
  }
  Q_gen(0, 0) = 0.25 * Q(0, 0) * dx(0) * dy(0);
  Q_gen(nx, 0) = 0.25 * Q(nx - 1, 0) * dx(nx - 1) * dy(0);
  Q_gen(0, ny) = 0.25 * Q(0, ny - 1) * dx(0) * dy(ny - 1);
  Q_gen(nx, ny) = 0.25 * Q(nx - 1, ny - 1) * dx(nx - 1) * dy(ny - 1);

  Array<double> T_old(nx + 1, ny + 1);

  int it = 0;
  double max_dif = 1e100;

  int max_it = 100000;
  double tol = 1e-9;
  double omega = 1.6;

  while (max_dif > tol && it < max_it) {
    it++;

    // Old temperature values
    T_old.fill(*T);

    (*T)(0, 0) = (0.5 * h * (dy(0) + dx(0)) * T_inf +
                  kdydx(0, 0) * (*T)(1, 0) +
                  kdxdy(0, 0) * (*T)(0, 1) +
                  Q_gen(0, 0)) / denom(0, 0);  // Bottom-left corner
    for (int i = 1; i < nx; i++)
      (*T)(i, 0) = (kdydx(i - 1, 0) * (*T)(i - 1, 0) +
                    kdydx(i, 0) * (*T)(i + 1, 0) +
                    h * dx(i) * T_inf +
                    kdxdy(i, 0) * (*T)(i, 1) +
                    Q_gen(i, 0)) / denom(i, 0);  // Bottom boundary
    (*T)(nx, 0) = (0.5 * (h * (dy(0) + dx(nx - 1)) * T_inf) +
                   kdydx(nx - 1, 0) * (*T)(nx - 1, 0) +
                   kdxdy(nx, 0) * (*T)(nx, 1) +
                   Q_gen(nx, 0)) / denom(nx, 0);  // Bottom-right corner
    for (int j = 1; j < ny; j++) {
      (*T)(0, j) = (h * dy(j) * T_inf +
                    kdydx(0, j) * (*T)(1, j) +
                    kdxdy(0, j - 1) * (*T)(0, j - 1) +
                    kdxdy(0, j) * (*T)(0, j + 1) +
                    Q_gen(0, j)) / denom(0, j);  // Right boundary
      for (int i = 1; i < nx; i++)
        (*T)(i, j) = (kdydx(i - 1, j) * (*T)(i - 1, j) +
                      kdydx(i, j) * (*T)(i + 1, j) +
                      kdxdy(i, j - 1) * (*T)(i, j - 1) +
                      kdxdy(i, j) * (*T)(i, j + 1) +
                      Q_gen(i, j)) / denom(i, j);  // Internal nodes
      (*T)(nx, j) = (kdydx(nx - 1, j) * (*T)(nx - 1, j) +
                     h * dy(j) * T_inf +
                     kdxdy(nx, j - 1) * (*T)(nx, j - 1) +
                     kdxdy(nx, j) * (*T)(nx, j + 1) +
                     Q_gen(nx, j)) / denom(nx, j);  // Right boundary
    }
    (*T)(0, ny) = (0.5 * (h * (dy(ny - 1) + dx(0)) * T_inf) +
                   kdydx(0, ny) * (*T)(0 + 1, ny) +
                   kdxdy(0, ny - 1) * (*T)(0, ny - 1) +
                   Q_gen(0, ny)) / denom(0, ny);  // Top-left corner
    for (int i = 1; i < nx; i++)
      (*T)(i, ny) = (kdydx(i - 1, ny) * (*T)(i - 1, ny) +
                     kdydx(i, ny) * (*T)(i + 1, ny) +
                     kdxdy(i, ny - 1) * (*T)(i, ny - 1) +
                     h * dx(i) * T_inf +
                     Q_gen(i, ny)) / denom(i, ny);  // Top boundary
    (*T)(nx, ny) = (0.5 * (h * (dy(ny - 1) + dx(nx - 1)) * T_inf) +
                    kdydx(nx - 1, ny) * (*T)(nx - 1, ny) +
                    kdxdy(nx, ny - 1) * (*T)(nx, ny - 1) +
                    Q_gen(nx, ny)) / denom(nx, ny);  // Top-right corner

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

  std::cout << "Number of iterations: " << it <<
            ", difference: " << max_dif << std::endl;
}

int main(int argc, char** argv) {
  // Parameters
  double Lx = 2;  // Length in x-direction [m]
  double Ly = 1;  // Length in x-direction [m]
  int nx = 64;    // Number of spatial regions in x-direction
  int ny = 32;    // Number of spatial regions in y-direction

  int n = nx * ny;      // Total number of spatial regions
  double dx = Lx / nx;  // x-direction interval
  double dy = Ly / ny;  // y-direction interval
  double A = dx * dy;   // Area of spatial region

  // x-direction nodes [m]
  Array<double> x(nx + 1, 1);
  for (int i = 0; i < nx + 1; i++)
    x(i) = i * dx;
  //x.print(4, 0);

  // y-direction nodes [m]
  Array<double> y(1, ny + 1);
  for (int j = 0; j < ny + 1; j++)
    y(j) = j * dy;
  //y.print(4, 0);

  // Thermal conductivity [W/m-K]
  double k0 = 1000;
  double k1 = 100000;
  Array<double> k(nx, ny);
  k.fill(k0);
  for (int i = nx / 2; i < nx; i++)
    for (int j = 3 * ny / 8; j < 5 * ny / 8; j++)
      k(i, j) = k1;
  //k.print(4, 0);

  // Volumetric heat source [W/m^3]
  double Q0_fwd = 1000000;  // Linear heat source [W/m]
  Array<double> Q_fwd(nx, ny);
  //Q_fwd(nx * 3 / 4, ny * 3 / 4) = Q0_fwd / A;
  for (int j = 0; j < ny; j++)
    Q_fwd(0, j) = Q0_fwd / A;
  //Q_fwd.print(4, 0);

  // Volumetric adjoint source [-]
  double Q0_adj = 1;  // Linear adjoint source [m^2]
  Array<double> Q_adj(nx, ny);
  //Q_adj(nx / 4, ny / 4) = Q0_adj / A;
  Q_adj(nx - 1, ny / 2) = Q0_adj;
  Q_adj(nx - 1, ny / 2 + 1) = Q0_adj;
  //Q_adj.print(4, 0);

  // Solve
  Array<double> T_fwd(nx + 1, ny + 1);
  Array<double> T_adj(nx + 1, ny + 1);
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  solve_fourier_2D(x, y, k, Q_fwd, 10000, 0, &T_fwd);
  solve_fourier_2D(x, y, k, Q_adj, 10000, 0, &T_adj);
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Elapsed time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
            << " ms" << std::endl;

  // Write forward and adjoing temperature to file
  T_fwd.printsci(4, "T_fwd.txt");
  T_adj.printsci(4, "T_adj.txt");

  return 0;
}

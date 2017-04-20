#include "Array.hpp"

void solve_fourier_2D(Array<double> x, Array<double> y, Array<double> k,
                      Array<double> Q) {

  int nx = x.get_nx() - 1;
  int ny = y.get_ny() - 1;

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
  for (int i = 1; i < nx; i++)
    for (int j = 1; j < ny; j++)
      denom(i, j) = kdydx(i - 1, j) + kdydx(i, j) +
                    kdxdy(i, j - 1) + kdxdy(i, j);

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
  Q_gen(0, ny) = 0.25 * Q(0, ny - 1) * dx(0) * dy(ny - 1);
  Q_gen(nx, 0) = 0.25 * Q(nx - 1, 0) * dx(nx - 1) * dy(0);
  Q_gen(nx, ny) = 0.25 * Q(nx - 1, ny - 1) * dx(nx - 1) * dy(ny - 1);

  Array<double> T(nx + 1, ny + 1);
  Array<double> T_old(nx + 1, ny + 1);

  int it = 0;
  double max_dif = 1e100;

  double tol = 1e-6;
  double omega = 1.6;

  while (max_dif > tol) {
    it++;

    // Old temperature values
    T_old.fill(T);

    // Internal nodes
    for (int j = 1; j < ny; j++)
      for (int i = 1; i < nx; i++)
        T(i, j) = (kdydx(i - 1, j) * T(i - 1, j) + kdydx(i, j) * T(i + 1, j) +
                   kdxdy(i, j - 1) * T(i, j - 1) + kdxdy(i, j) * T(i, j + 1) +
                   Q_gen(i, j)) / denom(i, j);

    // Successive over-relaxation
    T *= omega;
    T += T_old * (1 - omega);

    // Convergence criterion
    double dif = 0;
    max_dif = 0;
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        dif = T(i, j) - T_old(i, j);
        if (dif > max_dif)
          max_dif = dif;
      }
    }
  }

  std::cout << "Number of iterations: " << it << std::endl;
  T.print(3, 0);
}

int main(int argc, char** argv) {
  // Parameters
  double Lx = 1;  // Length in x-direction [m]
  double Ly = 1;  // Length in x-direction [m]
  int nx = 25;    // Number of spatial regions in x-direction
  int ny = 25;    // Number of spatial regions in y-direction

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
  double k0 = 100;
  Array<double> k(nx, ny);
  k.fill(k0);

  // Volumetric heat source [W/m^3]
  double Q0_fwd = 100000;  // Linear heat source [W/m]
  Array<double> Q_fwd(nx, ny);
  Q_fwd(nx * 3 / 4, ny / 4) =  Q0_fwd / A;

  // Volumetric adjoint source [-]
  double Q0_adj = 1;  // Linear adjoint source [m^2]
  Array<double> Q_adj(nx, ny);
  Q_adj(nx / 4, ny * 3 / 4) = Q0_adj / A;

  solve_fourier_2D(x, y, k, Q_fwd);

  return 0;
}

#include <iostream>

#include "Array.hpp"

void solve_fourier_2D(Array* xx, Array* yy, Array *k, Array *Q) {
  Array *dx, *dy, *kdx, *kdy, *kdydx, *kdxdy, *denom, *Q_gen, *T, *T_old;

  int nx = xx->get_nx() - 1;
  int ny = yy->get_ny() - 1;

  // dx: x-direction mesh; 1D; on intervals
  dx = new Array(nx, 1);
  for (int i = 0; i < nx; i++)
    dx->setVal(i, 0, xx->getVal(i + 1, 0) - xx->getVal(i, 0));

  // dy: y-direction mesh; 1D; on intervals
  dy = new Array(1, ny);
  for (int j = 0; j < ny; j++)
    dy->setVal(0, j, yy->getVal(0, j + 1) - yy->getVal(0, j));

  // kdy: k times y-direction mesh: 2D; on dx intervals
  kdy = new Array(nx, ny + 1);
  for (int j = 1; j < ny; j++)
    for (int i = 0; i < nx; i++)  // Internal
      kdy->setVal(i, j, 0.5*(k->getVal(i, j - 1)*dy->getVal(j - 1) +
                             k->getVal(i, j)*dy->getVal(j)));

  for (int i = 0; i < nx; i++) {  // Bottom, top
    kdy->setVal(i, 0, 0.5*k->getVal(i, 0)*dy->getVal(0));
    kdy->setVal(i, ny, 0.5*k->getVal(i, ny - 1)*dy->getVal(ny - 1));
  }

  // kdx: k times x-direction mesh: 2D; on dy intervals
  kdx = new Array(nx + 1, ny);
  for (int j = 0; j < ny; j++)
    for (int i = 1; i < nx; i++)  // Internal
      kdx->setVal(i, j, 0.5*(k->getVal(i - 1, j)*dx->getVal(i - 1) +
                             k->getVal(i, j)*dx->getVal(i)));
  for (int j = 0; j < ny; j++) {  // Left, right
    kdx->setVal(0, j, 0.5*k->getVal(0, j)*dx->getVal(0));
    kdx->setVal(nx, j, 0.5*k->getVal(nx - 1, j)*dx->getVal(nx - 1));
  }

  // kdy/dx
  kdydx = new Array(nx, ny + 1);
  for (int j = 0; j < ny + 1; j++)
    for (int i = 0; i < nx; i++)
      kdydx->setVal(i, j, kdy->getVal(i, j)/dx->getVal(i));
  delete kdy, dx;

  // kdx/dy
  kdxdy = new Array(nx + 1, ny);
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx + 1; i++)
      kdxdy->setVal(i, j, kdx->getVal(i, j)/dy->getVal(j));
  delete kdx, dy;

  // Denominator
  denom = new Array(nx + 1, ny + 1);
  for (int i = 1; i < nx; i++)
    for (int j = 1; j < ny; j++)
      denom->setVal(i, j, kdydx->getVal(i - 1, j) + kdydx->getVal(i, j) +
                          kdxdy->getVal(i, j - 1) + kdxdy->getVal(i, j));

  // Heat source
  Q_gen = new Array(nx + 1, ny + 1);
  for (int j = 1; j < ny; j++) {
    for (int i = 1; i < nx; i++) {
      Q_gen->setVal(i, j, 0.25*(
          Q->getVal(i - 1, j - 1)*dx->getVal(i - 1)*dy->getVal(j - 1) +
          Q->getVal(i - 1, j)*dx->getVal(i - 1)*dy->getVal(j) +
          Q->getVal(i, j - 1)*dx->getVal(i)*dy->getVal(j - 1) +
          Q->getVal(i, j)*dx->getVal(i)*dy->getVal(j)));
    }
    Q_gen->setVal(0, j, 0.25*(
        Q->getVal(0, j - 1)*dx->getVal(0)*dy->getVal(j - 1) +
        Q->getVal(0, j)*dx->getVal(0)*dy->getVal(j)));
    Q_gen->setVal(nx, j, 0.25*(
        Q->getVal(nx - 1, j - 1)*dx->getVal(nx - 1)*dy->getVal(j - 1) +
        Q->getVal(nx - 1, j)*dx->getVal(nx - 1)*dy->getVal(j)));
  }
  for (int i = 1; i < nx; i++) {
    Q_gen->setVal(i, 0, 0.25*(
        Q->getVal(i - 1, 0)*dx->getVal(i - 1)*dy->getVal(0) +
        Q->getVal(i, 0)*dx->getVal(i)*dy->getVal(0)));
    Q_gen->setVal(i, ny, 0.25*(
        Q->getVal(i - 1, ny - 1)*dx->getVal(i - 1)*dy->getVal(ny - 1) +
        Q->getVal(i, ny - 1)*dx->getVal(i)*dy->getVal(ny - 1)));
  }
  Q_gen->setVal(0, 0,
      0.25*Q->getVal(0, 0)*dx->getVal(0)*dy->getVal(0));
  Q_gen->setVal(0, ny,
      0.25*Q->getVal(0, ny - 1)*dx->getVal(0)*dy->getVal(ny - 1));
  Q_gen->setVal(nx, 0,
      0.25*Q->getVal(nx - 1, 0)*dx->getVal(nx - 1)*dy->getVal(0));
  Q_gen->setVal(nx, ny,
      0.25*Q->getVal(nx - 1, ny - 1)*dx->getVal(nx - 1)*dy->getVal(ny - 1));
  //Q_gen->print(5, 0);

  T = new Array(nx + 1, ny + 1);
  T_old = new Array(nx + 1, ny + 1);

  int it = 0;
  double max_dif = 1e100;

  double tol = 1e-6;
  double omega = 1.6;

  while (max_dif > tol) {
    it++;

    // Old temperature values
    T_old->fill(T);

    // Internal nodes
    for (int j = 1; j < ny; j++)
      for (int i = 1; i < nx; i++)
        T->setVal(i, j, (kdydx->getVal(i - 1, j)*T->getVal(i - 1, j) +
                         kdydx->getVal(i, j)*T->getVal(i + 1, j) +
                         kdxdy->getVal(i, j - 1)*T->getVal(i, j - 1) +
                         kdxdy->getVal(i, j)*T->getVal(i, j + 1) +
                         Q_gen->getVal(i, j))/denom->getVal(i, j));

    // Successive over-relaxation
    *T *= omega;
    *T += (*T_old)*(1 - omega);

    // Convergence criterion
    double dif = 0;
    max_dif = 0;
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        dif = T->getVal(i, j) - T_old->getVal(i, j);
        if (dif > max_dif) max_dif = dif;
      }
    }
  }

  // Cleanup
  delete kdydx, kdxdy, denom, Q_gen, T_old;

  std::cout << "Number of iterations: " << it << std::endl;
  T->print(3, 0);
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

  Array *xx, *yy, *k, *Q_fwd, *Q_adj;

  // x-direction nodes [m]
  xx = new Array(nx + 1, 1);
  for (int i = 0; i < nx + 1; i++)
    xx->setVal(i, 0, i * dx);

  // y-direction nodes [m]
  yy = new Array(1, ny + 1);
  for (int j = 0; j < ny + 1; j++)
    yy->setVal(0, j, j * dy);

  // Thermal conductivity [W/m-K]
  double k0 = 100;
  k = new Array(nx, ny);
  k->fill(k0);

  // Volumetric heat source [W/m^3]
  double Q0_fwd = 100000;  // Linear heat source [W/m]
  Q_fwd = new Array(nx, ny);
  Q_fwd->setVal(nx * 3 / 4, ny / 4, Q0_fwd / A);

  // Volumetric adjoint source [-]
  double Q0_adj = 1;  // Linear adjoint source [m^2]
  Q_adj = new Array(nx, ny);
  Q_adj->setVal(nx / 4, ny * 3 / 4, Q0_adj / A);

  solve_fourier_2D(xx, yy, k, Q_fwd);

  // Cleanup
  delete k, Q_fwd, Q_adj;

  return 0;
}

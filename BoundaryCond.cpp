#include "BoundaryCond.hpp"
#include "Common.hpp"

BoundaryCond::BoundaryCond() {
  type = BC_FIXED_TEMP;
  T_fixed = 0.;
  q_fixed = 0.;
  T_inf = 0.;
  h = 0.;
}

void BoundaryCond::set_type(int t, double p1, double p2) {
  type = t;
  if (t == BC_FIXED_TEMP)
    T_fixed = p1;
  else if (t == BC_FIXED_FLUX)
    q_fixed = p1;
  else if (t == BC_CONVECTIVE) {
    h = p1;
    T_inf = p2;
  }
  else if (t == BC_ADIABATIC) {}
  else
    throw_error("Unsupported boundary condition");
}

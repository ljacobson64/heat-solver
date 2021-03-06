#ifndef BOUNDARYCOND_H
#define BOUNDARYCOND_H

#include "Common.hpp"

#define BC_LEFT 0
#define BC_RIGHT 1
#define BC_BOTTOM 2
#define BC_TOP 3

#define BC_SIDE_MIN 0
#define BC_SIDE_MAX 3

#define BC_FIXED_TEMP 0
#define BC_FIXED_FLUX 1
#define BC_CONVECTIVE 2
#define BC_ADIABATIC 3

#define BC_TYPE_MIN 0
#define BC_TYPE_MAX 3

class BoundaryCond {
 public:
  BoundaryCond();  // Constructor
  ~BoundaryCond() {};  // Destructor

  int get_type() const {return type;};
  double get_T_fixed() const {return T_fixed;};
  double get_q_fixed() const {return q_fixed;};
  double get_T_inf() const {return T_inf;};
  double get_h() const {return h;};

  void set_type(int, double = 0., double = 0.);

 private:
  int type;
  double T_fixed;  // Fixed temperature
  double q_fixed;  // Fixed heat flux
  double T_inf;  // Ambient temperature [K]
  double h;  // Heat transfer coefficient [W/m^2-K]
};

typedef std::vector<BoundaryCond> BoundaryConds;

#endif

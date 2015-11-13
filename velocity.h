/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_VELOCITY_H
#define PS_VELOCITY_H

#include "pointers.h"

namespace PDPS_NS {

class Velocity : protected Pointers {
 public:

  Velocity(class PDPS *);
  void command(int, char **);
  void create(double, int);

 private:
  int dist_flag;
  int gid;
  int groupbit;
  class Compute *temperature;
  void set(int, char **);
  void scale(int, char **);
  void rescale(double, double);
};

}

#endif

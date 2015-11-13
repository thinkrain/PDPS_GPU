/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_RUN_H
#define PS_RUN_H

#include "pointers.h"

namespace PDPS_NS {

class Run : protected Pointers {
 public:
  Run(class PDPS *);
  void command(int, char **);
};

}

#endif

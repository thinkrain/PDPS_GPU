/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_OUTPUT_H
#define PS_OUTPUT_H

#include "pointers.h"

namespace PDPS_NS {

class Output : protected Pointers {
 public:
  int ndumps;
  class Dump **dump;

  // thermo
  class Thermo *thermo;        // Thermodynamic computations
  int thermo_every;            // output freq for thermo, 0 if first/last only

  Output(class PDPS *);
  ~Output();
  void init();
  void setup();

  void add_dump(int, char **);
  void modify_dump(int, char **);
  void delete_dump(char *id);

  void write();
  void print(const char *);

  void set_thermo(int, char **);     // set thermo output freqquency
  void create_thermo(int, char **);  // create a thermo style

private:
	int maxdump;
};

}

#endif

/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_CREATE_BOX_H
#define PS_CREATE_BOX_H

#include "pointers.h"

namespace PDPS_NS {

class CreateBox : protected Pointers {
 public:
	CreateBox(class PDPS *);
	void command(int, char **);
};

}

#endif

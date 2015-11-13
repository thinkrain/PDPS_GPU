/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_POST_PROCESSOR_H
#define PS_POST_PROCESSOR_H

#include "pointers.h"

namespace PDPS_NS {

class PostProcessor : protected Pointers {
public:

	PostProcessor(class PDPS *);
	void finalize();                    // provide information at the end of the simulation

private:
	int procid, nprocs;
	void analyze_particles();
	void analyze_neighbors();
	

};
}

#endif

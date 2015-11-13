/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(stress,ComputeStress)

#else

#ifndef PS_COMPUTE_STRESS_H
#define PS_COMPUTE_STRESS_H

#include "compute.h"
#include "fix_ave_spatial.h"

namespace PDPS_NS {

class ComputeStress : public Compute {
public:
	ComputeStress(class PDPS *, int, char **);
	virtual ~ComputeStress();
	void init();
	double compute_scalar();
	void compute_vector();

protected:
	int dim;
	int counter;
	char *style;
	int stress_flag;
	double boltz,nktv2p,volume;

	FILE *file;
	char *fname;

	int nevery;

	int fix_ave_id;
	FixAveSpatial *fixavespatial;
	int cid;
	double stress[6];
	double virial[6];

	void virial_compute(int);
	void Irving_Kirkwood();
	int check_bins(int);
};

}

#endif
#endif

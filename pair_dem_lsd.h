/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(dem/lsd,PairDEMLsd)

#else

#ifndef PS_PAIR_DEM_LSD_H
#define PS_PAIR_DEM_LSD_H

#include "pair.h"

namespace PDPS_NS {

class PairDEMLsd : public Pair {
public:

	PairDEMLsd(class PDPS *);
	~PairDEMLsd();

	void setup();
	void compute(int, int);
	void init_one(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);

protected:
	double cut_global;
	double **kn, **Cn;            // normal force
	double **kt, **Ct;            // tangential force
	double **mu;                  // friction coefficient
	double e;                     // restitution coefficient

	int tbsize;                   // page size for the pair_map in the pair_list, 
	                              // Could be the same as the one in the neighbor class

	void allocate();

private: 

	int rot_flag;                 // rotation flag

	void set_pairlist();
	void tags2str(char *, int, int);
};

}

#endif
#endif

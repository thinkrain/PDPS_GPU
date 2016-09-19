/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_CUDA_SPHPARTICLE_H
#define PS_CUDA_SPHPARTICLE_H

#include "pointers.h"
//#include "particle.h"
#include "cuda_engine.h"
#include "pdps_cuda.h"

namespace PDPS_NS {

class CUDASPHParticle : protected Pointers {
public:
	

	CUDASPHParticle(class PDPS *);
	~CUDASPHParticle();

	static double LoadFactor;
	int   HashTableSize;
	int  *devMapArray;
	unsigned int *devHashKey;
	unsigned int *devHashVal;
	unsigned int *devTry;

	virtual void SetMap();
	void data_atoms(int n, char *buf);

	
};

}

#endif


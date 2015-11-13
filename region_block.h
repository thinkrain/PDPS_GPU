/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS

RegionStyle(block, RegionBlock)

#else

#ifndef PS_REGION_BLOCK_H
#define PS_REGION_BLOCK_H

#include "region.h"

namespace PDPS_NS  {

class RegionBlock : public Region {
public:
	RegionBlock(class PDPS *, int, char **);
	~RegionBlock();

	int inside(double *);

	double find_distance(double *);

private:
	double xlo, xhi, ylo, yhi, zlo, zhi;
};

}

#endif
#endif

/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS

RegionStyle(point, RegionPoint)

#else

#ifndef PS_REGION_POINT_H
#define PS_REGION_POINT_H

#include "region.h"

namespace PDPS_NS  {

class RegionPoint : public Region {
public:

	RegionPoint(class PDPS *, int, char **);
	~RegionPoint();

	int inside(double *);                    // for point: only check whether particle coincide with some point

	double find_distance(double *);

	int find_point(char *);
	
private:
	
	
};

}

#endif
#endif

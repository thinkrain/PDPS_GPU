/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS

RegionStyle(intersect, RegionIntersect)

#else

#ifndef PS_REGION_INTERSECT_H
#define PS_REGION_INTERSECT_H

#include "region.h"

namespace PDPS_NS  {

	class RegionIntersect : public Region {
	public:
		RegionIntersect(class PDPS *, int, char **);
		~RegionIntersect();

		int inside(double *);

		double find_distance(double *);

	private:
		int *rid_list;         // region id list
		int nregions;          // # of regions to be stored


	};

}

#endif
#endif

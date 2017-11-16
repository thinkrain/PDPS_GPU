/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS

RegionStyle(polygon, RegionPolygon)

#else

#ifndef PS_REGION_POLYGON_H
#define PS_REGION_POLYGON_H

#include "region.h"

namespace PDPS_NS  {

class RegionPolygon : public Region {
public:
	RegionPolygon(class PDPS *, int, char **);
	~RegionPolygon();

	class RegionLine **edges;                         // edges

	void dynamic_check();

	int projection_inside(double *);                  // whether the point's projection is inside the plane
	int projection_inside(double *, double);          // whether the sphere's projection is inside the plane 
	int inside(double *);                             // wether a point is in this plane

	double find_distance(double *);                       // distance from a point to the plane
	double find_interaction_distance(double *, double *); // distance from a sphere to the plane
	void find_vector_to_rot_axis(double *, double *);

	void plane3D_map2_plane2D(double **);
	


private:
	FILE *test;           // debug purpose

	int region_id;

	double inv_abc;                                   // 1.0 / sqrt(a*a + b*b + c*c)

	int nps_initial;

	void calc_normal();
	double calc_area();
	void find_extent_bound();
	double find_projection_point(double *, double *, double *);
};

}

#endif
#endif

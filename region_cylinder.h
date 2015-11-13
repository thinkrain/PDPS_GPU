/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS

RegionStyle(cylinder, RegionCylinder)

#else

#ifndef PS_REGION_CYLINDER_H
#define PS_REGION_CYLINDER_H

#include "region.h"

namespace PDPS_NS  {

class RegionCylinder : public Region {
public:
	RegionCylinder(class PDPS *, int, char **);
	~RegionCylinder();

	int inside(double *);

	void dynamic_check();

	double find_distance(double *);                               // find distance between point to cylinder's nearest side surface
	double find_interaction_distance(double *, double *);         
	double find_projection_point(double *, double *, double *);   // find point's projection on cylinder's side surface
	double find_distance_to_rot_axis(double *);
	void find_vec_to_rot_axis(double *, double *);

private:

	double origin1[3], origin2[3];
	double a1, b1, c1, d1;                         // plane for the lower surface (origin1[3])
	double a2, b2, c2, d2;                         // plane for the lower surface (origin2[3])

	
	class RegionLine *axis;                        // axis, created by RegionLine class

	double find_distance_lower(double *);
	double find_distance_upper(double *);

	void find_extent_bound();
	void calc_normal();
	void rotate_cylinder(double **);
};

}

#endif
#endif

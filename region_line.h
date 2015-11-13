/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS

RegionStyle(line, RegionLine)

#else

#ifndef PS_REGION_LINE_H
#define PS_REGION_LINE_H

#include "region.h"

namespace PDPS_NS  {

class RegionLine : public Region {
public:

	RegionLine(class PDPS *, int, char **);
	~RegionLine();

	void dynamic_check();
	void calc_normal();

	int projection_inside(double *);                               // whether the point's projection is on the line
	int inside(double *);                                          // for line: only check whether particles are on the same line

	double find_distance(double *);                                // find distance from a point to the line
	double find_interaction_distance(double *, double *);  // find distance from a sphere to the line
	void find_vector(double *, double *);                          // find vector from a point to the line

	void rotate_point_around_axis(double *, double);               // rotate a point around a 3D axis   

	double find_projection_point(double *, double *, double *);    // find its projection point on the line, and return distance
	
private:

	//FILE *test;
	double calc_length();
	void find_extent_bound();
};

}

#endif
#endif

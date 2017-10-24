/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_REGION_H
#define PS_REGION_H

#include "pointers.h"

namespace PDPS_NS {

class Region : protected Pointers {
public:
	char *name, *style;

	int point_flag, line_flag, plane_flag, volume_flag;
	int rotate_flag;
	int dynamic_flag;

	double extent_lo[3], extent_hi[3], extent_le[3];
	double extent_xlo,extent_xhi;                        // bounding box on region
	double extent_ylo,extent_yhi;
	double extent_zlo,extent_zhi;
	double extent_xle, extent_yle, extent_zle;          

    int nps;                                             // No. of points
	double **coords;                                     // coordinates of the points
	double **v_coords;                                   // velocity of coordinate

	// translation
	// sin:  acc = A*sin(wt*t + phi)
	double init_time, cur_time;
	
	double tra_A[3];
	double tra_wt[3];
	double tra_phi[3];

	// rotation
	double romega, omega[3];
	int start, stable_start, stable_end, stop;     // rotation start, stable, and stop step
	int start_flag, stop_flag;
	double domega, omega_target;                   // increment of the rotation speed
	//double v[3];
	double naxis[3];
	class Region *rot_axis;                               // rotation axis, created by RegionLine class

	Region(class PDPS *, int, char **);
	virtual ~Region();

	void init();

	virtual void dynamic_check() {}

	virtual int inside(double *) = 0;
	virtual int inside(double *, double) { return 1; }
	int inside_a_plane(double *, int, double **);       // wether a point is in a plane
	
	virtual double find_distance(double *) = 0;
	virtual double find_interaction_distance(double *, double *) { return 0.0; }
	virtual void find_vector(double *, double *) {}
	double find_distance_to_a_point(double *, double *);
	double find_distance_to_a_line(double *, double *, double *);

	virtual double find_projection_point(double *, double *, double *) { return 0.0; }
	void find_projection_point_on_a_line(double *, double *, double *, double *, double *);

	virtual void rotate_point_around_axis(double *, double) {}

	//------------- point -----------------
	char **p_name;                                                    // point name

	virtual int find_point(char *) {return 0;}

	//------------- line -----------------
	double length;

	//------------- plain -----------------
	
	// polygon
	double a, b, c, d;                                // ax + by + cz = d
	double area;
	double **rot_coords;                                              // rotated coordinates after mapping the 3D plane into 2D plane (x-y)            
	double rot[3][3], inv_rot[3][3];                                  // rotation matrix and its inverse
	double normal[3];

	virtual int projection_inside(double *) {return 1;}
	virtual int projection_inside(double *, double) { return 1; }     // whether a sphere's projection is inside (include the edge case)
	virtual void plane3D_map2_plane2D(double **) {}
	
	//------------- volume -----------------
	double volume;

	// cylinder
	double radius, height; 

protected:
	int procid;

	void options(int, char **);
	void parse_point(char *, double *);

private:


};

}

#endif

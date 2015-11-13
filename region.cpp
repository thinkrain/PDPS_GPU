/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "error.h"
#include "parallel.h"
#include "phy_const.h"
#include "psmath.h"
#include "region.h"

using namespace PDPS_NS;
using namespace PhyConst;
using namespace PsMath_NS;

#define EPSILON 1.0e-6

enum{POINT,LINE,PLANE,VOLUME};

/* ---------------------------------------------------------------------- */

Region::Region(PDPS *ps, int narg, char **arg) : Pointers(ps)
{
	if (narg <= 2) error->all(FLERR, "Illegal region command");

	dynamic_flag = 0;
	start = stable_start = stable_end = stop = 0;
	start_flag = stop_flag = 0;

	// Region's name
	name = NULL;
	int n = strlen(arg[0]) + 1;
	name = new char[n];
	strcpy(name,arg[0]);

	// Region's style
	style = NULL;
	n = strlen(arg[1]) + 1;
	style = new char[n];
	strcpy(style,arg[1]);

	nps = 0;

	coords = NULL;
	v_coords = NULL;
	rot_coords = NULL;
	p_name = NULL;
	rot_axis = NULL;

	// initialize velocity to zero
	for (int i = 0; i < 3; i++) {
		//v[i] = 0.0;
		omega[i] = 0.0;
		romega = 0.0;
	}

	dynamic_flag = 0;
	rotate_flag = 0;

	procid = parallel->procid;
}

/* ---------------------------------------------------------------------- */

Region::~Region()
{
	delete [] name;
	delete [] style;
	name = NULL; 
    style = NULL;
}

/* ---------------------------------------------------------------------- */

void Region::init()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int Region::inside_a_plane(double *x0, int npoints, double **p)
{
	int i, j;
	int flag;

	double p1_norm, p2_norm, p1p2;
	double p1[3], p2[3];
	double c[3];
	double cos_theta, theta;
	theta = 0.0;

	for (i = 0; i < npoints; i++) {
		for (int k = 0; k < 3; k++) {
			p1[k] = p[i][k] - x0[k];	
		}

		j = (i+1) % nps;
		for (int k = 0; k < 3; k++) {
			p2[k] = p[j][k] - x0[k];	
		}

		p1_norm = Vec_Norm2(p1);
		p2_norm = Vec_Norm2(p2);
		p1p2 = p1_norm * p2_norm;
		if (p1p2 < EPSILON) {
			flag = 1;
			return flag;
		}
		else {
			cos_theta = Vec_Dot_Prod_3D(p1, p2) / p1p2;
			if (fabs(cos_theta) > 1.0) {
				if (fabs(cos_theta) - 1.0 < EPSILON) cos_theta = (cos_theta > 0.0 ? 1.0 : -1.0);
				else error->all(FLERR, "|cos_theta| is larger than 1.0, please debug the code");
			}
			Vec_Cross_Prod_3D(c, p1, p2);
			if (Vec_Dot_Prod_3D(c, normal) > 0) theta += acos(cos_theta);
			else theta -= acos(cos_theta);
		}
	}
	
	if (fabs(fabs(theta) - 2*PI) < EPSILON) flag = 1;
	else if (fabs(theta) < EPSILON) flag = 0;
	else {
		char str[128];
		sprintf(str, "The sum of the angle is equal to %f in radians and %f in degrees\n", theta, theta*180/PI);
		error->all(FLERR, str);
	}

  return flag;
}

/* ---------------------------------------------------------------------- */

void Region::parse_point(char *str, double *x0)
{
	int rid, pid;
	char r_name[128], p_name[128];

	char *ptr = strchr(str, '_');
	int n = ptr - str;

	// region name
	strncpy(r_name, str, n);
	r_name[n] = '\0';
	rid = domain->find_region(r_name);
	if (rid == -1) {
		char str1[128];
		sprintf(str1, "Cannot find the region id %s", r_name);
		error->all(FLERR, str1);
	}
	Region *region = domain->regions[rid];

	// point name
	ptr++;
	strcpy(p_name, ptr);
	
	pid = region->find_point(p_name);
	if (pid == -1) error->all(FLERR, "Cannot find the point id");

	x0[0] = region->coords[pid][0];
	x0[1] = region->coords[pid][1];
	x0[2] = region->coords[pid][2];
}

/* ----------------------------------------------------------------------
   Find the calculated distance direction
   For line, it is perpendicular as its normal directoin
   For plane, it is the same as its normal direction
   For volume, it depends on its shape
------------------------------------------------------------------------- */

void Region::find_projection_point_on_a_line(double *pro_coord, double *n, double *x0, double *p0, double *p1)
{
	double a[3], b[3];
	double dist;

	a[0] = x0[0] - p0[0];
	a[1] = x0[1] - p0[1];
	a[2] = x0[2] - p0[2];
	for (int i = 0; i < 3; i++) b[i] = p1[i] - p0[i];

	length = Vec_Norm2(b);

	double temp = Vec_Dot_Prod_3D(a, b);

	for (int i = 0; i < 3; i++) pro_coord[i] = temp / length * normal[i] + coords[0][i];

	for (int i = 0; i < 3; i++) n[i] = x0[i] - pro_coord[i];

	dist = Vec_Norm2(n);
	if (dist < EPSILON) {
		// the point is on the line
		for (int i = 0; i < 3; i++) n[i] = 0.0;
	}
	else {
		for (int i = 0; i < 3; i++) n[i] = n[i] / dist;
	}
}

/* ----------------------------------------------------------------------
   Find the distance between a point to a line;
   and the line is defined by one point and its direction
------------------------------------------------------------------------- */

double Region::find_distance_to_a_point(double *x0, double *p0)
{
	double r[3];
	for (int i = 0; i < 3; i++) r[i] = x0[i] - p0[i];

	double dist = Vec_Norm2(r);

	return dist; 
}

/* ----------------------------------------------------------------------
   Find the distance between a point to a line;
   and the line is defined by one point and its direction
------------------------------------------------------------------------- */

double Region::find_distance_to_a_line(double *x0, double *p0, double *n)
{
	double a[3], c[3];
	
	for (int i = 0; i < 3; i++) a[i] = x0[i] - p0[i];

	// c = a cross_prod b
	Vec_Cross_Prod_3D(c, a, n);
	
	double dist = Vec_Norm2(c);             // area2 = 2 * (area of the traingle)

	return dist;
}

/* ---------------------------------------------------------------------- */
/*
double Region::find_distance_to_rot_axis(double *x0)
{
	double dist = find_distance_to_a_line(x0, coords[0], normal);

	return dist;
}*/

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
#include "output.h"
#include "phy_const.h"
#include "psmath.h"
#include "region_line.h"
#include "update.h"

#include "particle.h"

using namespace PDPS_NS;
using namespace PhyConst;
using namespace PsMath_NS;

#define EPSILON 1.0e-10

enum{NONE, ROT, TRA_ACC_SIN};

/* ---------------------------------------------------------------------- */

RegionLine::RegionLine(PDPS *ps, int narg, char **arg) : Region(ps, narg, arg)
{
	if ((narg - 2) < 2) error->all(FLERR, "Illegal polygon line command");
	
	nps = 2;

	point_flag = 0;
	line_flag = 1;
	plane_flag = 0;
	volume_flag = 0;

	coords = new double *[nps];
	v_coords = new double *[nps];
	rot_coords = new double *[nps];
	for (int i = 0; i < nps; i++) {
		coords[i] = new double[3];
		v_coords[i] = new double[3];
		rot_coords[i] = new double[3];
	}

	// parse command
	int iarg = 2;
	int ip = 0;
	int rid, pid;
	// parse region point
	while (ip < 2) {
		parse_point(arg[iarg], coords[ip]);
		iarg++;
		ip++;
	}
	// parse options
	while (iarg < narg) {
		if (!strcmp(arg[iarg], "rot")) {
			rotate_flag = 1;
			dynamic_flag = ROT;
			int rid = domain->find_region(arg[iarg+1]);
			if (rid == -1) error->all(FLERR, "Illegal region id");
			rot_axis = domain->regions[rid];
			romega = atof(arg[iarg+2]);
			for (int i = 0; i < 3; i++) omega[i] = romega*rot_axis->normal[i];
			iarg += 3;
			while (iarg < narg) {
				if (!strcmp(arg[iarg], "start")) {
					start_flag = 1;
					start = atoi(arg[iarg+1]);
					stable_start = atoi(arg[iarg+2]);
					omega_target = romega;
					romega = 0.0;
					iarg += 3;
				}
				else if (!strcmp(arg[iarg], "stop")) {
					stop_flag = 1;
					stable_end = atoi(arg[iarg+1]);
					stop = atoi(arg[iarg+2]);
					iarg += 3;
				}
				else break;
			}
		}
		else if (!strcmp(arg[iarg], "tra_acc_sin")) {
			dynamic_flag = TRA_ACC_SIN;
			cur_time = 0.0;
			iarg++;
			for (int i = 0; i < 3; i++) {
				tra_A[i] = 0.0;
				tra_wt[i] = 0.0;
				tra_phi[i] = 0.0;
			}
			while (iarg < narg) {
				if (!strcmp(arg[iarg], "x")) {
					tra_A[0] = atof(arg[iarg+1]);
					tra_wt[0] = atof(arg[iarg+2]);
					tra_phi[0] = atof(arg[iarg+3]);
					iarg += 4;
				}
				else if (!strcmp(arg[iarg], "y")) {
					tra_A[1] = atof(arg[iarg+1]);
					tra_wt[1] = atof(arg[iarg+2]);
					tra_phi[1] = atof(arg[iarg+3]);
					iarg += 4;
				}
				else if (!strcmp(arg[iarg], "z")) {
					tra_A[2] = atof(arg[iarg+1]);
					tra_wt[2] = atof(arg[iarg+2]);
					tra_phi[2] = atof(arg[iarg+3]);
					iarg += 4;
				}
				else if (strcmp(arg[iarg], "start") == 0) {
					start_flag = 1;
					start = atoi(arg[iarg+1]);
					iarg += 2;
				}
				else if (strcmp(arg[iarg], "stop") == 0) {
					stop_flag = 1;
					stop = atoi(arg[iarg+1]);
					iarg += 2;
				}
				else break;
			}
		} // if (!strcmp(arg[iarg], "tra_acc_sin"))
		else error->all(FLERR, "Illegal region line options");
	}

	for (int i = 0; i < nps; i++)
	for (int j = 0; j < 3; j++) {
		v_coords[i][j] = 0.0;
	}

	length = calc_length();
	if (length < EPSILON) error->all(FLERR, "The line is too short, check your input again");

	calc_normal();
	find_extent_bound();

	//test = fopen("tran_analyze.txt", "w");
}

/* ---------------------------------------------------------------------- */

RegionLine::~RegionLine()
{
	for (int i = 0; i < nps; i++) {
		delete[] coords[i];
		coords[i] = NULL;
		delete[] rot_coords[i];
		rot_coords[i] = NULL;
	}
	delete[] coords;
	coords = NULL;
	delete[] rot_coords;
	rot_coords = NULL;
	delete[] v_coords;
	v_coords = NULL;
}

/* ----------------------------------------------------------------------
   Dynamic check:
   update rotation speed
   update rotation angle
------------------------------------------------------------------------- */

void RegionLine::dynamic_check()
{
	double cos_beta, sin_beta;
	double dist, dt;
	double theta;

	if (dynamic_flag == ROT) {
		// update rotating velocity
		int ntimestep = update->ntimestep;
		if (start_flag) {
			if (ntimestep < start) romega = 0.0;
			else if (ntimestep >= start && ntimestep < stable_start) {
				domega = omega_target / (stable_start - start);
				romega = domega * (ntimestep - start);
			}
			else if (ntimestep >= stable_start) {
				romega = omega_target;
			}
		}
		if (stop_flag) {
			if (ntimestep >= stop) {
				romega = 0.0;
				dynamic_flag = NONE;
			}
			else if (ntimestep >= stable_end && ntimestep < stop) {
				domega = omega_target / (stop - stable_end);
				romega = omega_target - domega * (ntimestep - stable_end);
			}
		}

		for (int i = 0; i < 3; i++) omega[i] = romega*rot_axis->normal[i];

		// rotate theta angle
		dt = update->dt;
		theta = romega * dt;
		for (int i = 0; i < nps; i++) {
			rot_axis->rotate_point_around_axis(coords[i], theta);
		}
	}
	else if (dynamic_flag == TRA_ACC_SIN) {
		double acc[3];
		double dt;

		if (stop_flag) {
			if (update->ntimestep >= stop) {
				dynamic_flag = NONE;
			}
		}
		if (update->ntimestep >= start && ((stop_flag == 1 && update->ntimestep <= stop) || stop_flag == 0)) {
			dt = update->dt;
			for (int i = 0; i < 3; i++) {
				acc[i] = tra_A[i] * sin(tra_wt[i] * cur_time + tra_phi[i]);
			}

			// translate the coordinate of points
			for (int i = 0; i < nps; i++)
			for (int j = 0; j < 3; j++) {
				v_coords[i][j] += acc[j] * dt;
				coords[i][j] += v_coords[i][j] * dt;
			}

			// This part is used to visualize the translation.
			// I did not delete it is because we may need it to visualize it for debug use
			/*for (int i = 0; i < particle->nlocal; i++) {
				if (particle->type[i] == 2 && !strcmp(name, "bottom")) {
					for (int j = 0; j < 3; j++) {
						particle->v[i][j] += acc[j] * dt;
						particle->x[i][j] += particle->v[i][j] * dt;
					}
				}
			}*/

			cur_time += dt;
		}
	}

	calc_normal();
	find_extent_bound();
}

/* ---------------------------------------------------------------------- */

int RegionLine::projection_inside(double *x0)
{
	double pro_coord[3], n[3];

	find_projection_point(pro_coord, n, x0);

	int ans = inside(pro_coord);
	
	return ans;
}

/* ----------------------------------------------------------------------
   If the point is on the line: dist = 0 && (x0 - p0).(x0 - p1) < 0
   or, the point is not on the line
------------------------------------------------------------------------- */

int RegionLine::inside(double *x0)
{
	double dist = find_distance(x0);

	if (fabs(dist) > EPSILON) {
		return 0;
	}

	double dot, a[3], b[3];
	for (int i = 0; i < 3; i++) {
		a[i] = x0[i] - coords[0][i];
		b[i] = x0[i] - coords[1][i];
	}

	dot = Vec_Dot_Prod_3D(a, b);
	
	if (dot <= 0.0) return 1;
	else return 0;
}

/* ----------------------------------------------------------------------
   Calculate the distance from a point to a line
------------------------------------------------------------------------- */

double RegionLine::find_distance(double *x0)
{
	double a[3], b[3], c[3];
	
	for (int i = 0; i < 3; i++) {
		a[i] = x0[i] - coords[0][i];
		b[i] = x0[i] - coords[1][i];
	}

	// c = a cross_prod b
	Vec_Cross_Prod_3D(c, a, b);
	
	double area2;             // area2 = 2 * (area of the traingle)
	area2 = Vec_Norm2(c);

	double dist = area2 / length;

	return dist;
}

/* ----------------------------------------------------------------------
   Calculate the distance from a sphere to the line:
   If the projection of the sphere's center is on the line, distance is
   the vector from the center to the line;
   Otherwise, the distance vector is the shortest distance to the
   head or tail
------------------------------------------------------------------------- */

double RegionLine::find_interaction_distance(double *n, double *x0) 
{
	int flag;
	double dist, pro_coord[3];

	dist = find_projection_point(pro_coord, n, x0);
	flag = inside(pro_coord);

	if (flag == 0) {
		// check the edges first
		dist = find_distance_to_a_point(x0, coords[0]);
		if (dist > EPSILON) {
			for (int i = 0; i < 3; i++) n[i] = (x0[i] - coords[0][i]) / dist;
		}
		else {
			for (int i = 0; i < 3; i++) n[i] = 0.0;
		}
		double temp;
		temp = find_distance_to_a_point(x0, coords[1]);
		if (dist > temp) {
			dist = temp;
			if (dist > EPSILON) {
				for (int i = 0; i < 3; i++) n[i] = (x0[i] - coords[1][i]) / dist;
			}
			else {
				for (int i = 0; i < 3; i++) n[i] = 0.0;
			}
		}
	}

	return dist;
}

/* ----------------------------------------------------------------------
   find the vector from point's projection to itself
------------------------------------------------------------------------- */

void RegionLine::find_vector(double *vec, double *x0)
{
	double pro_coord[3], n[3];
	find_projection_point(pro_coord, n, x0);
	for (int i = 0; i < 3; i++) vec[i] = x0[i] - pro_coord[i];
}

/* ----------------------------------------------------------------------
   find a point's projection vector:  x->x0
              x0
             / | \
         a  /  |  \
           /   |   \
    coords[0]--x----->coords[1]
	           b = coords[1] - coords[0]
	1. a.b = |a||b|cos(theta)
	2. |coords[0]->x| =|a|cos(theta)
	3. x = coords[0] + |coords[0]->x|*normal
------------------------------------------------------------------------- */

double RegionLine::find_projection_point(double *pro_coord, double *n, double *x0)
{
	double dist, a[3], b[3];

	for (int i = 0; i < 3; i++) {
		a[i] = x0[i] - coords[0][i];
		b[i] = coords[1][i] - coords[0][i];
	}

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
	return dist;
}

/* ---------------------------------------------------------------------- */

void RegionLine::calc_normal()
{
	for (int i = 0; i < 3; i++) {
		normal[i] = (coords[1][i]-coords[0][i]) / length;
	}
}

/* ---------------------------------------------------------------------- */

double RegionLine::calc_length()
{
	double r;

	r = 0;

	for (int i = 0; i < 3; i++) {
		r += (coords[1][i]-coords[0][i]) * (coords[1][i]-coords[0][i]);
	}

	r = sqrt(r);

	return r;
}

/* ---------------------------------------------------------------------- */

void RegionLine::find_extent_bound()
{
	// outbound of the region
	for (int i = 0; i < 3; i++) {
		extent_lo[i] = coords[0][i];
		extent_hi[i] = coords[0][i];
	}

	for (int i = 1; i < nps; i++) {
		for (int j = 0; j < 3; j++) {
			extent_lo[j] = MIN(extent_lo[j], coords[i][j]);
			extent_hi[j] = MAX(extent_hi[j], coords[i][j]);
		}
	}

	for (int i = 0; i < 3; i++) {
		extent_le[i] = extent_hi[i] - extent_lo[i];
	}

	extent_xlo = extent_lo[0];
	extent_xhi = extent_hi[0];
	extent_xle = extent_le[0];
	extent_ylo = extent_lo[1];
	extent_yhi = extent_hi[1];
	extent_yle = extent_le[1];
	extent_zlo = extent_lo[2];
	extent_zhi = extent_hi[2];
	extent_zle = extent_le[2];
}

/* ----------------------------------------------------------------------
   See the document "Rotation about an arbitrary axis in 3D dimensions"
------------------------------------------------------------------------- */

void RegionLine::rotate_point_around_axis(double *x0, double theta)
{
	double a, b, c;
	double u, v, w;
	double x, y, z;
	double cost, sint;
	double A[3];

	a = coords[0][0];
	b = coords[0][1];
	c = coords[0][2];
	u = normal[0];
	v = normal[1];
	w = normal[2];

	x = x0[0];
	y = x0[1];
	z = x0[2];

	cost = cos(theta);
	sint = sin(theta);
	x0[0] = (a*(v*v + w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1 - cost) +
		    x*cost + (-c*v + b*w - w*y + v*z)*sint;
	x0[1] = (b*(u*u + w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1 - cost) +
		    y*cost + (c*u - a*w + w*x - u*z)*sint;
	x0[2] = (c*(u*u + v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1 - cost) +
		    z*cost + (-b*u + a*v - v*x + u*y)*sint;
}

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
#include "particle.h"
#include "phy_const.h"
#include "psmath.h"
#include "region_cylinder.h"
#include "region_line.h"
#include "update.h"

using namespace PDPS_NS;
using namespace PhyConst;
using namespace PsMath_NS;

#define EPSILON 1e-10

enum{NONE, ROT, TRA_ACC_SIN};

/* ---------------------------------------------------------------------- */

RegionCylinder::RegionCylinder(PDPS *ps, int narg, char **arg) : Region(ps, narg, arg)
{
	point_flag = 0;
	line_flag = 0;
	plane_flag = 0;
	volume_flag = 1;

	dynamic_flag = NONE;
	nps = 2;

	axis = NULL;

	coords = new double *[nps];
	v_coords = new double *[nps];
	rot_coords = new double *[nps];
	for (int i = 0; i < nps; i++) {
		coords[i] = new double[3];
		v_coords[i] = new double[3];
		rot_coords[i] = new double[3];
	}

	// parse points
	int iarg = 2;
	int ip = 0;
	while (ip < 2) {
		parse_point(arg[iarg], coords[ip]);
		iarg++;
		ip++;
	}
	radius = atof(arg[iarg++]);

	height = (coords[1][0]-coords[0][0])*(coords[1][0]-coords[0][0]) +
		     (coords[1][1]-coords[0][1])*(coords[1][1]-coords[0][1]) +
		     (coords[1][2]-coords[0][2])*(coords[1][2]-coords[0][2]);
	height = sqrt(height);
	if (height < EPSILON) error->all(FLERR, "The cylinder's height is too small");

	volume = PI * radius * radius * height;

	// parse options
	while (iarg < narg) {
		if (!strcmp(arg[iarg], "rotate")) {
			rotate_flag = ROT;
			int rid = domain->find_region(arg[iarg+1]);
			if (rid == -1) error->all(FLERR, "Illegal region id");
			rot_axis = domain->regions[rid];
			omega_target = atof(arg[iarg+2]);
			romega = omega_target;
			iarg += 3;
			while (iarg < narg) {
				if (strcmp(arg[iarg], "start") == 0) {
					start_flag = 1;
					start = atoi(arg[iarg+1]);
					stable_start = atoi(arg[iarg+2]);
					iarg += 3;
				}
				else if (strcmp(arg[iarg], "stop") == 0) {
					stop_flag = 1;
					stable_end = atoi(arg[iarg+1]);
					stop = atoi(arg[iarg+2]);
					iarg += 3;
				}
				else break;
			}
			for (int i = 0; i < 3; i++) omega[i] = omega_target*rot_axis->normal[i];
		} // if (!strcmp(arg[iarg], "rotate"))
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
		} // else if (!strcmp(arg[iarg], "tra_acc_sin"))
		else error->all(FLERR, "Illegal region cylinder options");
	}

	// create cylinder axis
	int narg_line = narg - 1;
	char **arg_line = new char*[narg_line];
	for (int i = 0; i < narg_line; i++) arg_line[i] = new char[128];
	sprintf(arg_line[0], "%s-cylinder-axis", arg[0]);
	sprintf(arg_line[1], "line");
	strcpy(arg_line[2], arg[2]);
	strcpy(arg_line[3], arg[3]);
	for (int i = 4; i < narg_line; i++) {
		strcpy(arg_line[i], arg[i+1]);
	}
	axis = new RegionLine(ps, narg_line, arg_line);
	for (int i = 0; i < narg_line; i++) delete[] arg_line[i];
	delete[] arg_line;

	// initialize translational and rotation velocity
	for (int i = 0; i < nps; i++)
	for (int j = 0; j < 3; j++) {
		v_coords[i][j] = 0.0;
	}

	if (dynamic_flag == ROT) {
		if (start_flag) {
			if (start > stable_start) error->all(FLERR, "Illegal region clylinder rotating start options");
		}
		if (stop_flag) {
			if (stable_end > stop) error->all(FLERR, "Illegal region clylinder rotating stop options");
		}
	}

	calc_normal();
	rotate_cylinder(rot_coords);
	find_extent_bound();
}

/* ---------------------------------------------------------------------- */

RegionCylinder::~RegionCylinder()
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

	delete axis;
	axis = NULL;
}

/* ---------------------------------------------------------------------- */

void RegionCylinder::dynamic_check()
{
	double dt, theta;
	
	if (dynamic_flag == ROT) {
		int ntimestep = update->ntimestep;

		if (start_flag == 0) romega = omega_target;
		else {
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

		if (update->ntimestep >= start && ((stop_flag == 1 && update->ntimestep < stop) || stop_flag == 0)) {
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
			for (int i = 0; i < particle->nlocal; i++) {
				if (particle->type[i] == 3 && !strcmp(name, "vessel")) {
					for (int j = 0; j < 3; j++) {
						particle->v[i][j] += acc[j] * dt;
						particle->x[i][j] += particle->v[i][j] * dt;
					}
				}
			}
			cur_time += dt;
		}
	}

	// update the axis
	axis->dynamic_check();

	calc_normal();
	find_extent_bound();
}

/* ---------------------------------------------------------------------- */

void RegionCylinder::calc_normal()
{
	normal[0] = (coords[1][0] - coords[0][0]) / height;
	normal[1] = (coords[1][1] - coords[0][1]) / height;
	normal[2] = (coords[1][2] - coords[0][2]) / height;

	a1 = normal[0];
	b1 = normal[1];
	c1 = normal[2];
	d1 = a1*coords[0][0] + b1*coords[0][1] + c1*coords[0][2];

	a2 = -normal[0];
	b2 = -normal[1];
	c2 = -normal[2];
	d2 = a2*coords[1][0] + b2*coords[1][1] + c2*coords[1][2];
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside the cylinder bound
   inside = 0 if x,y,z is outside the cylinder bound
------------------------------------------------------------------------- */

int RegionCylinder::inside(double *x0)
{
	double dist = find_distance(x0);

	if (dist < 0) return 0;

	double dist1 = find_distance_lower(x0);

	if (dist1 >= 0 && dist1 <= height) return 1;
	else return 0;
}

/* ----------------------------------------------------------------------
   calculate distance to the surface: 
   >= 0: closer than the surface to the axis
   < 0: further than the surface to the axis
------------------------------------------------------------------------- */

double RegionCylinder::find_distance(double *x0)
{
	double a[3], b[3], c[3];
	
	for (int i = 0; i < 3; i++) {
		a[i] = x0[i] - coords[0][i];
		b[i] = x0[i] - coords[1][i];
	}

	// c = a cross_prod b
	Vec_Cross_Prod_3D(c, a, b);
	double area2;

	area2 = Vec_Norm2(c);       // area2 = 2 * (area of the traingle)

	double dist = radius - area2 / height;

	return dist;
}

/* ----------------------------------------------------------------------
   Find interaction distance from a sphere to the cylinder
   The edge case has not been developed yet
------------------------------------------------------------------------- */

double RegionCylinder::find_interaction_distance(double *n, double *x0)
{
	int flag;
	double dist, pro_coord[3];

	dist = find_projection_point(pro_coord, n, x0);
	flag = inside(pro_coord);
	// finish it in the future
	if (flag == 0) {
		/*............................*/
		/*............................*/
		/*............................*/
	}

	return dist;
}

/* ---------------------------------------------------------------------- */

double RegionCylinder::find_distance_lower(double *x0)
{
	double dist = a1*x0[0] + b1*x0[1] + c1*x0[2] - d1;

	return dist;
}

/* ----------------------------------------------------------------------
   Find four points coordinate for the lower and upper surface
------------------------------------------------------------------------- */

void RegionCylinder::find_extent_bound()
{
	double pro_p1[4][3], pro_p2[4][3];
	double p1[4][3], p2[4][3];
	double shiftx, shifty;

	shiftx = radius;
	shifty = radius;

	pro_p1[0][0] = rot_coords[0][0] + shiftx;
	pro_p1[0][1] = rot_coords[0][1] + shifty;
	pro_p1[0][2] = rot_coords[0][2]; 

	pro_p1[1][0] = rot_coords[0][0] - shiftx;
	pro_p1[1][1] = rot_coords[0][1] + shifty;
	pro_p1[1][2] = rot_coords[0][2];

	pro_p1[2][0] = rot_coords[0][0] - shiftx;
	pro_p1[2][1] = rot_coords[0][1] - shifty;
	pro_p1[2][2] = rot_coords[0][2];
	
	pro_p1[3][0] = rot_coords[0][0] + shiftx;
	pro_p1[3][1] = rot_coords[0][1] - shifty;
	pro_p1[3][2] = rot_coords[0][2];

	pro_p2[0][0] = rot_coords[1][0] + shiftx;
	pro_p2[0][1] = rot_coords[1][1] + shifty;
	pro_p2[0][2] = rot_coords[1][2]; 

	pro_p2[1][0] = rot_coords[1][0] - shiftx;
	pro_p2[1][1] = rot_coords[1][1] + shifty;
	pro_p2[1][2] = rot_coords[1][2];

	pro_p2[2][0] = rot_coords[1][0] - shiftx;
	pro_p2[2][1] = rot_coords[1][1] - shifty;
	pro_p2[2][2] = rot_coords[1][2];
	
	pro_p2[3][0] = rot_coords[1][0] + shiftx;
	pro_p2[3][1] = rot_coords[1][1] - shifty;
	pro_p2[3][2] = rot_coords[1][2];

	for (int i = 0; i < 4; i++) {
		Matrix_Prod_3D(p1[i], inv_rot, pro_p1[i]);
		Matrix_Prod_3D(p2[i], inv_rot, pro_p2[i]);
	}

	extent_lo[0] = extent_hi[0] = p1[0][0];
	extent_lo[1] = extent_hi[1] = p1[0][1];
	extent_lo[2] = extent_hi[2] = p1[0][2];
	
	for (int i = 1; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			extent_lo[j] = MIN(extent_lo[j], p1[i][j]);
			extent_hi[j] = MAX(extent_hi[j], p1[i][j]);
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			extent_lo[j] = MIN(extent_lo[j], p2[i][j]);
			extent_hi[j] = MAX(extent_hi[j], p2[i][j]);
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
   Find the projection point on the cylinder side surface
   The direction returned is from the projected point to the point x0
------------------------------------------------------------------------- */

double RegionCylinder::find_projection_point(double *pro_coord, double *n, double *x0)
{
	double pro_coord1[3];
	double r[3];
	double dist;

	dist = axis->find_projection_point(pro_coord1, n, x0);
	if (dist < EPSILON) {
		return radius;
	}

	for (int i = 0; i < 3; i++) r[i] = radius * n[i];
	for (int i = 0; i < 3; i++) pro_coord[i] = pro_coord1[i] + r[i];
	for (int i = 0; i < 3; i++) {
		n[i] = x0[i] - pro_coord[i];
	}
	dist = Vec_Norm2(n);
	if (dist < EPSILON) {
		for (int i = 0; i < 3; i++) n[i] = 0.0;
	}
	else {
		for (int i = 0; i < 3; i++) n[i] = n[i] / dist;
	}

	return dist;
}

/* ----------------------------------------------------------------------
   Rotate axis to project a 3D plane into a 2D plane.
   The new coordinate system is presented by x', y', z', where its origin
   is fixed, z' is aligned with the normal vector n of the 3D plane 
   (note that vector x', y', and z' are still presented in old coordinate 
   system), and x' is chosen along the line by connecting the first and 
   second point stored in the coords[][].
   Hence, we will have the rotation matrix R as: 
   |1 0 0|         |x'[0] y'[0] z'[0]|
   |0 1 0| = rot * |x'[1] y'[1] z'[1]|
   |0 0 1|         |x'[2] y'[2] z'[2]|
   or B = rot * A
   It means after the rotation, the vector x', y', and z' in the new 
   coordinate system, will be unit vector (1 0 0)', (0 1 0)', and (0 0 1)
   Therefore, one can use R to find all new coordinates of points (coords[][]) 
   under the x'-y'-z' system. Sunch values are stored in p_x[][]
------------------------------------------------------------------------- */

void RegionCylinder::rotate_cylinder(double **p_coords)
{
	int i, j;
	double A[3][3], B[3][3], inv_A[3][3];
	double vec0[3], vec1[3], vec2[3];

	// ----------- Construct matrix A ---------------

	// align z  to normal direction of the plane 
	A[0][2] = normal[0];
	A[1][2] = normal[1];
	A[2][2] = normal[2];

	// align x to the line connectin the first two points
	double temp1, temp2[3], inv_temp1;
	double p_assist[3], tx, tz;

	// (tx, 0, tz).(normal[0], normal[1], normal[2]) = 0
	// tx^2 + tz^2 = 1.0
	// tx^2(normal[0]^2 + normal[2]^2) = normal[2]*normal[2]

	double sq_xz = (normal[0]*normal[0] + normal[2]*normal[2]);
	if (fabs(sq_xz) < EPSILON) {
		tx = 1.0;
		tz = 0.0;
	}
	else {
		tx = normal[2]*normal[2] / sq_xz;
		tz = sqrt(1.0 - tx);
		tx = sqrt(tx);
	}

	p_assist[0] = coords[0][0] + tx*radius;
	p_assist[1] = coords[0][1];
	p_assist[2] = coords[0][2] + tz*radius;

	temp1 = 0.0;
	for (i = 0; i < 3; i++) {
		temp1 += (p_assist[i] - coords[0][i]) * (p_assist[i] - coords[0][i]);
	}
	temp1 = sqrt(temp1);
	inv_temp1 = 1.0 / temp1;

	for (i = 0; i < 3; i++) temp2[i] = (p_assist[i] - coords[0][i]) * inv_temp1;

	for (i = 0; i < 3; i++) A[i][0] = temp2[i];

	// y = x cross_prod z

	for (i = 0; i < 3; i++) {
		vec0[i] = A[i][0]; 
		vec2[i] = A[i][2];
	}

	Vec_Cross_Prod_3D(vec1, vec0, vec2);
	// y'
	for (i = 0; i < 3; i++) A[i][1] = vec1[i];

	// x'
	for (i = 0; i < 3; i++) A[i][0] = vec0[i];
	// z'
	for (i = 0; i < 3; i++) A[i][2] = vec2[i];

	for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++) {
		inv_rot[i][j] = A[i][j];
	}

	// ----------- Construct matrix B ---------------

	for (i = 0; i < 3; i++) 
	for (j = 0; j < 3; j++) {
		B[i][j] = 0.0;
	}
	B[0][0] = 1.0;
	B[1][1] = 1.0;
	B[2][2] = 1.0;

	// ---------- calculate rotation matrix R -------------

	int success = Matrix_Inverse_3D(inv_A, A);
	if (success == 0) error->all(FLERR, "The determinant of the matrix is close to 0");

	Matrix_Prod_3D(rot, B, inv_A);

	// --------- project all points to the new coordinate system x'-y'-z' -----------

	double vec[3];
	double p_vec[3];

	for (i = 0; i < nps; i++) {
		for (j = 0; j < 3; j++) vec[j] = coords[i][j];

		Matrix_Prod_3D(p_vec, rot, vec);

		for (j = 0; j < 3; j++) p_coords[i][j] = p_vec[j];
	}
}

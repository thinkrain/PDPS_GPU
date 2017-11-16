/* ----------------------------------------------------------------------
    < Particle Dynamics Parallel Simulator (PDPS) >
	Copyright(C) <2014>  <Author: Lingqi Yang>
	Email: ly2282@columbia.edu

	This program is free software : you can redistribute it and / or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include <iostream>
#include "domain.h"
#include "error.h"
#include "force.h"
#include "fix_wall_dem.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"
#include "pair_dem_lsd.h"
#include "parallel.h"
#include "particle.h" 
#include "pair_list.h"
#include "psmath.h"
#include "region.h"
#include "update.h"
#include "modify.h"

#include "output.h"

#include "pdps_cuda.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "device_functions.h"

using namespace PDPS_NS;
using namespace FixConst;
using namespace PsMath_NS;
#define EPSILON 1.0e-10
#define PI 3.1416
#define float_double double 
#define Max_Wall 12

enum{LSD};
enum{BOUND, REGION};

/* ---------------------------------------------------------------------- */

FixWallDEM::FixWallDEM(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 12) error->all(FLERR,"Illegal fix wall/dem command");

	nwalls = 0;

	wall_rid = NULL;
	wall_flag = NULL;


	if (!strcmp(arg[3], "lsd")) {
		wall_style = LSD;
	}
	else error->all(FLERR, "Illegal fix wall/dem style");

	// parse options
	int iarg = 4;
	
	if ((narg-7-4) % 2 != 0) error->all(FLERR, "Illegal fix wall/dem option");
	int nwalls_initial = (narg-7-4) / 2;

	wall_rid = new int[nwalls_initial];
	wall_flag = new int[nwalls_initial];
	//vel_wall = new double*[nwalls_initial];
	/*for (int i = 0; i < nwalls_initial; i++) {
		vel_wall[i] = new double[3];
		for (int j = 0; j < 3; j++) vel_wall[i][j] = 0.0;
	}*/
	
	while (iarg < narg) {
		if (!strcmp(arg[iarg], "xlo")) {
			wall[nwalls] = 0;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "xhi")) {
			wall[nwalls] = 1;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "ylo")) {
			wall[nwalls] = 2;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "yhi")) {
			wall[nwalls] = 3;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		} 
		else if (!strcmp(arg[iarg], "zlo")) {
			if (domain->dim == 2) error->all(FLERR, "It is a 2D simulation");
			wall[nwalls] = 4;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "zhi")) {
			if (domain->dim == 2) error->all(FLERR, "It is a 2D simulation");
			wall[nwalls] = 5;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "region")) {
			wall_flag[nwalls] = REGION;
			wall_rid[nwalls] = domain->find_region(arg[iarg+1]);
			if (wall_rid[nwalls] == -1) error->all(FLERR, "Cannot find the region id");
			nwalls++;
		}
		else break;
		iarg += 2;
	} // while (iarg < narg)


	if (narg - iarg != 7) error->all(FLERR, "Illegal fix wall/dem command");

	rot_flag = 1;
	e = atof(arg[iarg++]);
	kn = atof(arg[iarg++]);
	Cn = atof(arg[iarg++]);
	kt = atof(arg[iarg++]);
	Ct = atof(arg[iarg++]);
	mu = atof(arg[iarg++]);
	cut = atof(arg[iarg++]);

	rneigh = cut + neighbor->rskin;

	pair_list = NULL;
	tbsize = 10000;
	pair_list = new PairList(ps, tbsize);
	pair_list->ifields = 1;
	pair_list->dfields = 3;
	
	// assign GPU parameters
	//modify->gpunwalls = nwalls;
	//modify->gpuwall = wall;
	//modify->gpucoords0 = coords0;
	//modify->gpuwall_rid = wall_rid;
	//modify->gpuwall_flag = wall_flag;
	//modify->gpuwall_cut = cut;
	hostWallContactFlag = NULL;
	hostWallDrijtx = NULL;
	hostWallDrijty = NULL;
	hostWallDrijtz = NULL;
	devWallContactFlag = NULL;
	devWallDrijtx = NULL;
	devWallDrijty = NULL;
	devWallDrijtz = NULL;
	devWall = NULL;
	devCoords0 = NULL;
	devWall_rid = NULL;
	devWall_flag = NULL;
	//devCoord1X = NULL;
	//devCoord1Y = NULL;
	//devCoord1Z = NULL;
	//devCoord2X = NULL;
	//devCoord2Y = NULL;
	//devCoord2Z = NULL;
	//devCoord3X = NULL;
	//devCoord3Y = NULL;
	//devCoord3Z = NULL;
	//devCoord4X = NULL;
	//devCoord4Y = NULL;
	//devCoord4Z = NULL;
	//devStyle = NULL;
	//devRadiusCylinder = NULL;
	//devHeight = NULL;
	//devA = NULL;
	//devB = NULL;
	//devC = NULL;
	//devD = NULL;
}

/* ---------------------------------------------------------------------- */

FixWallDEM::~FixWallDEM()
{
	delete[] wall_rid;
	wall_rid = NULL;
	
	delete[] wall_flag;
	wall_flag = NULL;

	delete pair_list;
	pair_list = NULL;

	free(hostWallContactFlag);
	free(hostWallDrijtx);
	free(hostWallDrijty);
	free(hostWallDrijtz);
	cudaFree(devWallContactFlag);
	cudaFree(devWallDrijtx);
	cudaFree(devWallDrijty);
	cudaFree(devWallDrijtz);
	cudaFree(devWall);
	cudaFree(devCoords0);
	cudaFree(devWall_rid);
	cudaFree(devWall_flag);
	/*cudaFree(devCoord1X);
	cudaFree(devCoord1Y);
	cudaFree(devCoord1Z);
	cudaFree(devCoord2X);
	cudaFree(devCoord2Y);
	cudaFree(devCoord2Z);
	cudaFree(devCoord3X);
	cudaFree(devCoord3Y);
	cudaFree(devCoord3Z);
	cudaFree(devCoord4X);
	cudaFree(devCoord4Y);
	cudaFree(devCoord4Z);
	cudaFree(devStyle);
	cudaFree(devRadiusCylinder);
	cudaFree(devHeight);
	cudaFree(devA);
	cudaFree(devB);
	cudaFree(devC);
	cudaFree(devD);*/

}

/* ---------------------------------------------------------------------- */

void FixWallDEM::init()
{
	int size = particle->nlocal;
	pair_list->init_hash(size);

	// Initialize GPU parameters
	hostWallContactFlag = (int *)malloc(particle->nmax * nwalls * sizeof(int));
	hostWallDrijtx = (double *)malloc(particle->nmax * nwalls * sizeof(double));
	hostWallDrijty = (double *)malloc(particle->nmax * nwalls * sizeof(double));
	hostWallDrijtz = (double *)malloc(particle->nmax * nwalls * sizeof(double));

	cudaMalloc(&devWallContactFlag, particle->nmax * nwalls * sizeof(int));
	cudaMemset(devWallContactFlag, 0, particle->nmax * nwalls * sizeof(int));
	cudaMalloc(&devWallDrijtx, particle->nmax * nwalls * sizeof(double));
	cudaMalloc(&devWallDrijty, particle->nmax * nwalls * sizeof(double));
	cudaMalloc(&devWallDrijtz, particle->nmax * nwalls * sizeof(double));
	cudaMalloc(&devWall, Max_Wall * sizeof(int));
	cudaMalloc(&devCoords0, Max_Wall * sizeof(double));
	cudaMalloc(&devWall_rid, nwalls * sizeof(int));
	cudaMalloc(&devWall_flag, nwalls * sizeof(int));
	cudaMemcpy(devWall, wall, Max_Wall * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoords0, coords0, Max_Wall * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devWall_rid, wall_rid, nwalls * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devWall_flag, wall_flag, nwalls * sizeof(int), cudaMemcpyHostToDevice);

	//devStyle = NULL;
	//devRadiusCylinder = NULL;
	//devHeight = NULL;
	//devA = NULL;
	//devB = NULL;
	//devC = NULL;
	//devD = NULL;

	/*std::string keyCylinder = "cylinder";
	std::string keyPolygon = "polygon";
	for (int iwall = 0; iwall < nwalls; iwall++){
		int rid = wall_rid[iwall];
		if (!strncmp(>regions[i]->style, keyCylinder.c_str(), keyCylinder.size())){
			devStyle[iwall] = 1;
			devRadiusCylinder[iwall] = regions[rid]->radius;
			devHeight[iwall] = regions[rid]->height;
			devCoord1X[iwall] = regions[rid]->coords[0][0];
			devCoord1Y[iwall] = regions[rid]->coords[0][1];
			devCoord1Z[iwall] = regions[rid]->coords[0][2];
			devCoord2X[iwall] = regions[rid]->coords[1][0];
			devCoord2Y[iwall] = regions[rid]->coords[1][1];
			devCoord2Z[iwall] = regions[rid]->coords[1][2];
			devA[iwall] = regions[rid]->normal[0];
			devB[iwall] = regions[rid]->normal[1];
			devC[iwall] = regions[rid]->normal[2];
	}*/

}

/* ---------------------------------------------------------------------- */

void FixWallDEM::setup()
{
	set_pairlist();
}

/* ---------------------------------------------------------------------- */

int FixWallDEM::setmask()
{
	int mask = 0;
	mask |= PRE_FORCE;
	return mask;
}


__device__ void gpu_pre_force_dem_lsd(
    double *devCoordX, double *devCoordY, double *devCoordZ,
    double *devVeloX, double *devVeloY, double *devVeloZ,  
	int *devWallContactFlag, double *devWallDrijtx, double *devWallDrijty, double *devWallDrijtz,
	int *devWall, double *devCoords0, int *devWall_rid, int *devWall_flag,
	double *devMass, double *devRadius, int *devType, int *devMask, const int nlocal, const int groupbit, const double dt, 
	const int radius_flag, const int torque_flag,
	const double kn, const double Cn, const double kt, const double Ct, const double mu, const double cut, const int nwalls,
	double *devForceX, double *devForceY, double *devForceZ,
	const int iparticle, const double delta, double *n, const int iwall){
	
	//int table, index;
	float_double drijn, drijnx, drijny, drijnz;                // normal displacement
	float_double drijtx, drijty, drijtz;                       // tangential displacement
	float_double vijx, vijy, vijz;                             // relative velocity: vij = vj - vi
	float_double vijn, vijnx, vijny, vijnz;                    // relative velocity along normal direction
	float_double vijt, vijt_inv, vijtx, vijty, vijtz;          // relative velocity along tangential direction
	float_double fijn, fijnx, fijny, fijnz;                    // normal force
	float_double fijt, fijtx, fijty, fijtz;                    // tangential force	
	float_double vel_wall[3];
	//double fijt_a[3], fijt, fijtx, fijty, fijtz;         // tangential force	
	float_double nx, ny, nz;                                   // unit vector along normal direction
	float_double tx, ty, tz;                                   // unit vector along tangential direction
	float_double temp;
	//double omegainij[3], omegajnij[3];                   // omega_i cross nij
	//double torqueij[3];                                  // torque
	//double Li[3], Lj[3];                                 // vector from the center of particle to the contact point 
    //
	//double *radius = particle->radius;
	//double **omega = particle->omega;
	//double **torque = particle->torque;
    //
    //
	//Region **regions = domain->regions;
    //
	nx = n[0];
	ny = n[1];
	nz = n[2];
   
	// Overlap
	drijn = delta;
	drijnx = drijn * nx;
	drijny = drijn * ny;
	drijnz = drijn * nz;
	
	// may need to consider wall translational velocity in the future
	for (int i = 0; i < 3; i++) vel_wall[i] = 0.0;
    
	if (devWall_flag[iwall] == 1) {
		int rid = devWall_rid[iwall];
		//if (regions[rid]->dynamic_flag) {
		//	for (int j = 0; j < 3; j++) vel_wall[j] = regions[rid]->v_coords[0][j];
		//}
	}
	// vij = vi - vj + (wi cross Ri) - (wj cross Rj) (relative velocity)
	vijx = devVeloX[iparticle] - vel_wall[0];
	vijy = devVeloY[iparticle] - vel_wall[1];
	vijz = devVeloZ[iparticle] - vel_wall[2];
    
	//if (particle->radius_flag) {
	//	for (int i = 0; i < 3; i++) Li[i] = (radius[iparticle] - drijn)*(-n[i]);
	//	Vec_Cross_Prod_3D(omegainij, omega[iparticle], Li);
	//	vijx += omegainij[0];
	//	vijy += omegainij[1];
	//	vijz += omegainij[2];
	//	if (wall_flag[iwall] == REGION) {
	//		int rid = wall_rid[iwall];
	//		if (regions[rid]->rotate_flag) {
	//			double pro_coord[3], nn[3];
	//			regions[rid]->find_projection_point(pro_coord, nn, x[iparticle]);
	//			regions[rid]->rot_axis->find_vector(Lj, pro_coord);
	//			Vec_Cross_Prod_3D(omegajnij, regions[rid]->omega, Lj);
	//			// the direction of the center of the region to the contact point is opposite of n
	//			vijx += -omegajnij[0];
	//			vijy += -omegajnij[1];
	//			vijz += -omegajnij[2];
	//		}
	//	}
	//}
    
	// |vijn| = vij . nij
	vijn = vijx * nx + vijy * ny + vijz * nz;
    // vijn = |vijn| . nij
	vijnx = vijn*nx;
	vijny = vijn*ny;
	vijnz = vijn*nz;
    
	// vijt = vij - (vij . nij) . nij
	vijtx = vijx - vijnx;
	vijty = vijy - vijny;
	vijtz = vijz - vijnz;
	vijt = sqrt(vijtx*vijtx + vijty*vijty + vijtz*vijtz);
	if (fabs(vijt) < EPSILON) vijt_inv = 0.0;
	else vijt_inv = 1.0/vijt;
    
	// tij = tij/|tij|
	tx = vijtx * vijt_inv;
	ty = vijty * vijt_inv;
	tz = vijtz * vijt_inv;
		
	// calculate drijt
	// first time to contact
	int contact_flag = devWallContactFlag[iparticle * nwalls + iwall];
	if (contact_flag == 0) {
		if (fabs(vijn) < EPSILON) temp = dt;
		else temp = MIN(fabs(drijn/vijn),dt);
		drijtx = vijtx * temp;
		drijty = vijty * temp;
		drijtz = vijtz * temp;
		devWallContactFlag[iparticle * nwalls + iwall] = 1;
	}
	// during the same contact
	else {
		drijtx = devWallDrijtx[iparticle * nwalls + iwall];
		drijty = devWallDrijty[iparticle * nwalls + iwall];
		drijtz = devWallDrijtz[iparticle * nwalls + iwall];
		
		// update the tang. disp. for the next time step
		drijtx = drijtx + vijtx*dt;
		drijty = drijty + vijty*dt;
		drijtz = drijtz + vijtz*dt;
    
		temp = drijtx*nx + drijty*ny + drijtz*nz;
    
		drijtx = drijtx - temp*nx;
		drijty = drijty - temp*ny;
		drijtz = drijtz - temp*nz;
	}
	devWallDrijtx[iparticle * nwalls + iwall] = drijtx;
	devWallDrijty[iparticle * nwalls + iwall] = drijty;
	devWallDrijtz[iparticle * nwalls + iwall] = drijtz;
	
	// fijn = kn*delta_rijn - Cn*vijn
	fijnx = kn*drijnx - Cn*vijnx;
	fijny = kn*drijny - Cn*vijny;
	fijnz = kn*drijnz - Cn*vijnz;
	fijn = sqrt(fijnx*fijnx + fijny*fijny + fijnz*fijnz);

	// tangential force
	// fijt = -kt*delta_rijt - Ct*vijt
	fijtx = -kt * drijtx - Ct * vijtx;
	fijty = -kt * drijty - Ct * vijty;
	fijtz = -kt * drijtz - Ct * vijtz;
	fijt = sqrt(fijtx*fijtx + fijty*fijty + fijtz*fijtz); 
	temp = fabs(mu*fijn);
	// |fijt| > mu*|fnij|
	if (fijt > temp) {
		fijtx = -temp * tx;
		fijty = -temp * ty;
		fijtz = -temp * tz;
	}
    //
	devForceX[iparticle] += (fijnx + fijtx);
	devForceY[iparticle] += (fijny + fijty);
	devForceZ[iparticle] += (fijnz + fijtz);
    
	//if (particle->torque_flag) {
	//	fijt_a[0] = fijtx;
	//	fijt_a[1] = fijty;
	//	fijt_a[2] = fijtz;
	//	Vec_Cross_Prod_3D(torqueij, Li, fijt_a);
    //
	//	torque[iparticle][0] += torqueij[0];
	//	torque[iparticle][1] += torqueij[1];
	//	torque[iparticle][2] += torqueij[2];
	//}
}

__device__ double gpu_regionLine_projection(
    double *n, double *pro_coord, double *x0, double *coord1, double *coord2, double *normal, double length){
	
	double a[3], b[3], dist;
	
	for (int i = 0; i < 3; i++) {
		a[i] = x0[i] - coord1[i];
		b[i] = coord2[i] - coord1[i];
	}
	double temp = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	for (int i = 0; i < 3; i++) pro_coord[i] = temp / length * normal[i] + coord1[i];
	for (int i = 0; i < 3; i++) n[i] = x0[i] - pro_coord[i];
	dist = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	if (dist < EPSILON) {
		for (int i = 0; i < 3; i++) n[i] = 0.0;
	}
	else {
		for (int i = 0; i < 3; i++) n[i] = n[i] / dist;
	}
	
	return dist;
}

__device__ double gpu_regionLine_distance(
    double *n, double *x0, double *coord1, double *coord2, double *normal, double length){
	
	int flag = 1;
	double dist, dist1, pro_coord[3];
	double a[3], b[3], c[3], r[3];

	// find projected point on the line
	dist = gpu_regionLine_projection(n, pro_coord, x0, coord1, coord2, normal, length);
	
	// check if the point is insode the line segment
	for (int i = 0; i < 3; i++) {
		a[i] = pro_coord[i] - coord1[i];
		b[i] = pro_coord[i] - coord2[i];
	}
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
	double area2 = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
	dist1 = area2 / length;
	if (fabs(dist1) > EPSILON) {
		flag = 0;
	}
	if (flag == 1) {
		double dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
		if (dot <= 0.0) flag = 1;
		else flag = 0;
	}

	// if the point is not inside, compute the distance from the end points
	if (flag == 0) {
		for (int i = 0; i < 3; i++) r[i] = x0[i] - coord1[i];
		dist = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
		if (dist > EPSILON) {
			for (int i = 0; i < 3; i++) n[i] = (x0[i] - coord1[i]) / dist;
		}
		else {
			for (int i = 0; i < 3; i++) n[i] = 0.0;
		}
		for (int i = 0; i < 3; i++) r[i] = x0[i] - coord2[i];
		double temp = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
		if (dist > temp) {
			dist = temp;
			if (dist > EPSILON) {
				for (int i = 0; i < 3; i++) n[i] = (x0[i] - coord2[i]) / dist;
			}
			else {
				for (int i = 0; i < 3; i++) n[i] = 0.0;
			}
		}
	}

	return dist;
}
		
__device__ double gpu_regionCylinder_distance(
    double *n, double x, double y, double z,
	double devRadiusCylinder, double devHeight,
	double devCoord1X, double devCoord1Y, double devCoord1Z,
	double devCoord2X, double devCoord2Y, double devCoord2Z,
	double devA, double devB, double devC){
	
	// initialization
	double dist, a[3], b[3], x0[3], coord1[3], coord2[3], normal[3];
	double r[3], pro_coord[3], pro_coord1[3];
	x0[0] = x;
	x0[1] = y;
	x0[2] = z;
	coord1[0] = devCoord1X;
	coord1[1] = devCoord1Y;
	coord1[2] = devCoord1Z;
	coord2[0] = devCoord2X;
	coord2[1] = devCoord2Y;
	coord2[2] = devCoord2Z;
	normal[0] = devA;
	normal[1] = devB;
	normal[2] = devC;

	// projection point on the axis
	dist = gpu_regionLine_projection(n, pro_coord1, x0, coord1, coord2, normal, devHeight);
	if (dist < EPSILON) {
		return devRadiusCylinder;
	}
	
	// compute distance from the cylinder
	for (int i = 0; i < 3; i++) r[i] = devRadiusCylinder * n[i];
	for (int i = 0; i < 3; i++) pro_coord[i] = pro_coord1[i] + r[i];
	for (int i = 0; i < 3; i++) {
		n[i] = x0[i] - pro_coord[i];
	}
	dist = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	if (dist < EPSILON) {
		for (int i = 0; i < 3; i++) n[i] = 0.0;
	}
	else {
		for (int i = 0; i < 3; i++) n[i] = n[i] / dist;
	}
	
	return dist;
}


__device__ double gpu_regionPolygon_distance(
    double *n, double x, double y, double z,
	double devCoord1X, double devCoord1Y, double devCoord1Z,
	double devCoord2X, double devCoord2Y, double devCoord2Z,
	double devCoord3X, double devCoord3Y, double devCoord3Z,
	double devCoord4X, double devCoord4Y, double devCoord4Z,
	double devA, double devB, double devC, double devD){
	
	double dist, a[3], b[3], x0[3], coords[4][3];
	double r[3], pro_coord[3], pro_coord1[3], distVec[3], normal[3], p1[3], p2[3], theta;
	int flag = 0;
	//
	x0[0] = x;
	x0[1] = y;
	x0[2] = z;
	coords[0][0] = devCoord1X;
	coords[0][1] = devCoord1Y;
	coords[0][2] = devCoord1Z;
	coords[1][0] = devCoord2X;
	coords[1][1] = devCoord2Y;
	coords[1][2] = devCoord2Z;
	coords[2][0] = devCoord3X;
	coords[2][1] = devCoord3Y;
	coords[2][2] = devCoord3Z;
	coords[3][0] = devCoord4X;
	coords[3][1] = devCoord4Y;
	coords[3][2] = devCoord4Z;
	normal[0] = devA;
	normal[1] = devB;
	normal[2] = devC;
	theta = 0.;
	
	// compute distance from the plane
	dist = devA*x0[0] + devB*x0[1] + devC*x0[2] - devD; 
	for (int i = 0; i < 3; i++) {
		n[i] = dist / fabs(dist) * normal[i];
		distVec[i] = dist*normal[i];
		pro_coord[i] = x0[i] - distVec[i];
	}
	
	// check if we are inside the polygon
	int nps = 4; // We are limited to 4-node polygon
	for (int i = 0; i < nps; i++) {
		for (int k = 0; k < 3; k++) p1[k] = coords[i][k] - pro_coord[k];
		int j = (i+1) % nps;
		for (int k = 0; k < 3; k++) p2[k] = coords[j][k] - pro_coord[k];
		double p1_norm = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);
		double p2_norm = sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);
		double p1p2 = p1_norm * p2_norm;
		if (p1p2 < EPSILON) {
			flag = 1;
			break;
		}
		else {
			double cos_theta = (p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]) / p1p2;
			if (fabs(cos_theta) > 1.0) {
				if (fabs(cos_theta) - 1.0 < EPSILON) cos_theta = (cos_theta > 0.0 ? 1.0 : -1.0);
			}
			double c[3];
			c[0] = p1[1]*p2[2] - p1[2]*p2[1];
			c[1] = p1[2]*p2[0] - p1[0]*p2[2];
			c[2] = p1[0]*p2[1] - p1[1]*p2[0];
			if ((c[0]*normal[0] + c[1]*normal[1] + c[2]*normal[2]) > 0) theta += acos(cos_theta);
			else theta -= acos(cos_theta);
		}
	}
	if (flag == 0) {
		if (fabs(fabs(theta) - 2*PI) < 1e-3) flag = 1;
		else if (fabs(theta) < 1e-3) flag = 0;
	}
	// find the distance from the edges
	if (flag == 0) {
		double n1[3], length, normalLine[3];
		for (int i = 0; i < 3; i++) r[i] = coords[0][i] - coords[nps-1][i];
		length = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
		for (int i = 0; i < 3; i++) normalLine[i] = r[i]/length;
		dist = gpu_regionLine_distance(n, x0, coords[nps-1], coords[0], normalLine, length);
		double temp;
		for (int ip = 0; ip < nps-1; ip++) {
			for (int i = 0; i < 3; i++) r[i] = coords[ip+1][i] - coords[ip][i];
			length = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			for (int i = 0; i < 3; i++) normalLine[i] = r[i]/length;
			temp = gpu_regionLine_distance(n1, x0, coords[ip], coords[ip+1], normalLine, length);
			if (fabs(dist) > fabs(temp)) {
				dist = temp;
				for (int j = 0; j < 3; j++) n[j] = n1[j];
			}
		}
	}
	
	return dist;
}

__device__ void gpu_test(
	double *devCoordX, double *devCoordY, double *devCoordZ,
	double *devVeloX, double *devVeloY, double *devVeloZ, 
	int *devWallContactFlag, double *devWallDrijtx, double *devWallDrijty, double *devWallDrijtz,
	int *devWall, double *devCoords0, int *devWall_rid, int *devWall_flag,
	double *devMass, double *devRadius, int *devType, int *devMask, const int nlocal, const int groupbit, const double dt,
	const int radius_flag, const int torque_flag,
	const double kn, const double cn, const double kt, const double ct, 
	const double mu, const double cut, const int nwalls){
	//devWallDrijty[0] = 1.6;
}

__device__ void gpu_test3(const double dt
	){

}

__device__ void gpu_test2(
	double *devCoordX, double *devCoordY, double *devCoordZ,
	double *devVeloX, double *devVeloY, double *devVeloZ,
	int *devWallContactFlag, double *devWallDrijtx, double *devWallDrijty, double *devWallDrijtz,
	int *devWall, double *devCoords0, int *devWall_rid, int *devWall_flag,
	double *devMass, double *devRadius, int *devType, int *devMask, const int nlocal, const int groupbit, const double dt,
	const int radius_flag, const int torque_flag,
	const double kn, const double Cn, const double kt, const double Ct, const double mu, const double cut, const int nwalls,
	double *devForceX, double *devForceY, double *devForceZ,
	const int iparticle, const double delta, double *n, const int iwall){

}
/* ---------------------------------------------------------------------- */

__global__ void gpuPreforce(
    double *devCoordX, double *devCoordY, double *devCoordZ,
    double *devVeloX, double *devVeloY, double *devVeloZ,  
	int *devWallContactFlag, double *devWallDrijtx, double *devWallDrijty, double *devWallDrijtz,
	int *devWall, double *devCoords0, int *devWall_rid, int *devWall_flag,
	double *devMass, double *devRadius, int *devType, int *devMask, const int nlocal, const int groupbit, const double dt, 
	const int radius_flag, const int torque_flag,
	const double kn, const double Cn, const double kt, const double Ct, const double mu, const double cut, const int nwalls,
	double *devForceX, double *devForceY, double *devForceZ,
	int *devStyle, double *devRadiusCylinder, double *devHeight,
	double *devCoord1X, double *devCoord1Y, double *devCoord1Z,
	double *devCoord2X, double *devCoord2Y, double *devCoord2Z,
	double *devCoord3X, double *devCoord3Y, double *devCoord3Z,
	double *devCoord4X, double *devCoord4Y, double *devCoord4Z,
	double *devA, double *devB, double *devC, double *devD, double ntimestep){	
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	//float_double wf, ix, iy, iz, rsq, rij, radius_cut, rijx, rijy, rijz, q, tmp, fi, fj, irho, jrho, rij_inv, delVdotDelR;
	//float_double ivx, ivy, ivz, fviscx, fviscy, fviscz, wfd, fpair, imass, jmass, nx, ny, nz;
	//float_double drijn, drijnx, drijny, drijnz, vijx, vijy, vijz, Li, Lj;
	//float_double omegainijx, omegainijy, omegainijz, omegajnijx, omegajnijy, omegajnijz;
	//float_double vijn, vijnx, vijny, vijnz, vijt, vijtx, vijty, vijtz, vijt_inv, tx, ty, tz;
	//float_double drijtx, drijty, drijtz, temp, torqueijx, torqueijy, torqueijz;
	//float_double fijn, fijnx, fijny, fijnz, fijt, fijtx, fijty, fijtz;
	//unsigned int j, jj, itype, jtype, jnum;
	//int contact_flag;

	for (i = i; i < nlocal; i += blockDim.x * gridDim.x){
		for (int iwall = 0; iwall < nwalls; iwall++) {
			int rid;
			int dim, side;
			double dist, n[3];
			if (devMask[i] & groupbit) {
				if (devWall_flag[iwall] == 0) {
					dim = devWall[iwall] / 2;
					side = devWall[iwall] % 2;
					n[0] = 0.;
					n[1] = 0.;
					n[2] = 0.;
					if (side == 0) {
						n[dim] = 1.0;
						if (dim==0) dist = devCoordX[i] - devCoords0[iwall];
						else if (dim==1) dist = devCoordY[i] - devCoords0[iwall];
						else if (dim==2) dist = devCoordZ[i] - devCoords0[iwall];
					}
					else {
						n[dim] = -1.0;
						if (dim==0) dist = devCoords0[iwall] - devCoordX[i];
						else if (dim==1) dist = devCoords0[iwall] - devCoordY[i];
						else if (dim==2) dist = devCoords0[iwall] - devCoordZ[i];
					}
				}
				else if (devWall_flag[iwall] == 1) {
					rid = devWall_rid[iwall];
					//dist = domain->regions[rid]->find_interaction_distance(n, x[i]);
					if (devStyle[rid]==1) dist = gpu_regionCylinder_distance(
													n, devCoordX[i], devCoordY[i], devCoordZ[i],
													devRadiusCylinder[rid], devHeight[rid],
													devCoord1X[rid], devCoord1Y[rid], devCoord1Z[rid],
													devCoord2X[rid], devCoord2Y[rid], devCoord2Z[rid],
													devA[rid], devB[rid], devC[rid]);
					else if (devStyle[rid]==2) dist = gpu_regionPolygon_distance(
														n, devCoordX[i], devCoordY[i], devCoordZ[i],
														devCoord1X[rid], devCoord1Y[rid], devCoord1Z[rid],
														devCoord2X[rid], devCoord2Y[rid], devCoord2Z[rid],
														devCoord3X[rid], devCoord3Y[rid], devCoord3Z[rid],
														devCoord4X[rid], devCoord4Y[rid], devCoord4Z[rid],
														devA[rid], devB[rid], devC[rid], devD[rid]);
				}
				

				// check if the particle is within the interaction range
				if (radius_flag == 0 && fabs(dist) < cut) {
					gpu_pre_force_dem_lsd(devCoordX, devCoordY, devCoordZ,
						devVeloX, devVeloY, devVeloZ,  
						devWallContactFlag, devWallDrijtx, devWallDrijty, devWallDrijtz,
						devWall, devCoords0, devWall_rid, devWall_flag,
						devMass, devRadius, devType, devMask, nlocal, groupbit, dt, 
						radius_flag, torque_flag,
						kn, Cn, kt, Ct, mu, cut, nwalls,
						devForceX, devForceY, devForceZ,
						i, cut-fabs(dist), n, iwall);
				}
				else if (radius_flag == 1 && fabs(dist) < devRadius[i]) {
					gpu_pre_force_dem_lsd(devCoordX, devCoordY, devCoordZ,
						devVeloX, devVeloY, devVeloZ,  
						devWallContactFlag, devWallDrijtx, devWallDrijty, devWallDrijtz,
						devWall, devCoords0, devWall_rid, devWall_flag,
						devMass, devRadius, devType, devMask, nlocal, groupbit, dt, 
						radius_flag, torque_flag,
						kn, Cn, kt, Ct, mu, cut, nwalls,
						devForceX, devForceY, devForceZ,
						i, devRadius[i]-fabs(dist), n, iwall);
				}
				else {
					devWallContactFlag[i * nwalls + iwall] = 0;
					devWallDrijtx[i * nwalls + iwall] = 0.;
					devWallDrijty[i * nwalls + iwall] = 0.;
					devWallDrijtz[i * nwalls + iwall] = 0.;
				}
			} // if (mask[i] & groupbit)
		} // for (iwall = 0; iwall < nwalls; iwall++)
	} // for (i, nlocal)
}

	
/* ---------------------------------------------------------------------- */

void FixWallDEM::pre_force()
{
	//int iwall, itag;

	//double **x = particle->x;
	//double **v = particle->v;
	//double **f = particle->f;
	//int *mask = particle->mask;
	//int *tag = particle->tag;
	int radius_flag = particle->radius_flag;
	//double *radius = particle->radius;
	int nlocal = particle->nlocal;
	double dt = update->dt;
	//std::cout<<"update->ntimestep="<<update->ntimestep<<std::endl;

	// char str[BYTES];
	// int ipair;
    // 
	// int nflag = neighbor->nflag;
    // 
	// if (nflag == 1) {
	// 	set_pairlist();
	// }
    // 
	// PairList::HashPair **hash_pair = pair_list->hash_pair;
    // 
	// int rid;
	// int dim, side;
	// double dist, n[3];
    // 
	// for (iwall = 0; iwall < nwalls; iwall++) {
	// 	for (int i = 0; i < nlocal; i++) {
	// 		if (mask[i] & groupbit) {
	// 			itag = tag[i];
	// 			tag2str(str, itag, iwall);
	// 			ipair = pair_list->find_hash(str);
	// 			if (wall_flag[iwall] == BOUND) {
	// 				dim = wall[iwall] / 2;
	// 				side = wall[iwall] % 2;
	// 				for (int j = 0; j < 3; j++) n[j] = 0.0;
	// 				if (side == 0) {
	// 					dist = x[i][dim] - coords0[iwall];
	// 					n[dim] = 1.0;
	// 				}
	// 				else {
	// 					dist = coords0[iwall] - x[i][dim];
	// 					n[dim] = -1.0;
	// 				}
	// 			}
	// 			else if (wall_flag[iwall] == REGION) {			
	// 				rid = wall_rid[iwall];
	// 				dist = domain->regions[rid]->find_interaction_distance(n, x[i]);
	// 			}
    // 
	// 			// check if the particle is within the interaction range
	// 			if (radius_flag == 0 && fabs(dist) < cut) {
	// 				if (ipair == -1) ipair = pair_list->insert_hash(str);
	// 				pre_force_dem_lsd(ipair, i, cut-fabs(dist), n, iwall);
	// 			}
	// 			else if (radius_flag == 1 && fabs(dist) < radius[i]) {
	// 				if (ipair == -1) ipair = pair_list->insert_hash(str);
	// 				pre_force_dem_lsd(ipair, i, radius[i]-fabs(dist), n, iwall);
	// 			}
	// 			else {
	// 				if (ipair > -1) {
	// 					pair_list->set_zero(ipair);
	// 				}
	// 			}
	// 		} // if (mask[i] & groupbit)
	// 	} // for (int i = 0; i < nlocal; i++)
	// } // for (iwall = 0; iwall < nwalls; iwall++)
	cudaError_t error_t;

	gpuPreforce << < int(nlocal + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >
	 (particle->devCoordX, particle->devCoordY, particle->devCoordZ,
	  particle->devVeloX, particle->devVeloY, particle->devVeloZ,
	  devWallContactFlag, devWallDrijtx, devWallDrijty, devWallDrijtz,
	  devWall, devCoords0, devWall_rid, devWall_flag,
	  particle->devMass, particle->devRadius, particle->devType,
	  particle->devMask, nlocal, groupbit, dt, radius_flag, particle->torque_flag, 
	  kn, Cn, kt, Ct, mu, cut, nwalls,
	  particle->devForceX, particle->devForceY, particle->devForceZ,
	  domain->devStyle, domain->devRadiusCylinder, domain->devHeight,
	  domain->devCoord1X, domain->devCoord1Y, domain->devCoord1Z,
	  domain->devCoord2X, domain->devCoord2Y, domain->devCoord2Z,
	  domain->devCoord3X, domain->devCoord3Y, domain->devCoord3Z,
	  domain->devCoord4X, domain->devCoord4Y, domain->devCoord4Z,
	  domain->devA, domain->devB, domain->devC, domain->devD, update->ntimestep);

	
	//error_t = cudaMemcpy(neighbor->hostForceX, particle->devForceX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(neighbor->hostForceY, particle->devForceY, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(neighbor->hostForceZ, particle->devForceZ, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//
	//error_t = cudaMemcpy(neighbor->hostVeloX, particle->devVeloX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(neighbor->hostVeloY, particle->devVeloY, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(neighbor->hostVeloZ, particle->devVeloZ, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);

	//error_t = cudaMemcpy(hostWallContactFlag, devWallContactFlag, particle->nmax * nwalls * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(hostWallDrijtx, devWallDrijtx, particle->nmax * nwalls * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(hostWallDrijty, devWallDrijty, particle->nmax * nwalls * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(hostWallDrijtz, devWallDrijtz, particle->nmax * nwalls * sizeof(double), cudaMemcpyDeviceToHost);

	//error_t = cudaMemcpy(neighbor->hostForceX, particle->devForceX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(neighbor->hostForceZ, particle->devForceZ, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);

	//error_t = cudaMemcpy(domain->hostStyle, domain->devStyle, domain->nregions * sizeof(int), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostRadiusCylinder, domain->devRadiusCylinder, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostHeight, domain->devHeight, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord1X, domain->devCoord1X, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord1Y, domain->devCoord1Y, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord1Z, domain->devCoord1Z, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord2X, domain->devCoord2X, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord2Y, domain->devCoord2Y, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord2Z, domain->devCoord2Z, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord3X, domain->devCoord3X, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord3Y, domain->devCoord3Y, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord3Z, domain->devCoord3Z, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord4X, domain->devCoord4X, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord4Y, domain->devCoord4Y, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostCoord4Z, domain->devCoord4Z, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostA, domain->devA, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostB, domain->devB, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostC, domain->devC, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(domain->hostD, domain->devD, domain->nregions * sizeof(double), cudaMemcpyDeviceToHost);

	//if (hostWallDrijtx[0] < 0.5)
	//	error_t = error_t;
 //	 particle->TransferG2C();
}

/* ----------------------------------------------------------------------
   only reset pairlist when nflag = 1
------------------------------------------------------------------------- */

void FixWallDEM::set_pairlist()
{
	int ipair;
	int itag;
	int dim, side;
	int iwall;
	char str[BYTES];

	int nlocal = particle->nlocal;
	int radius_flag = particle->radius_flag;
	int *tag = particle->tag;
	int *mask = particle->mask;
	double **x = particle->x;
	double *radius = particle->radius;

	// initialize # of undetermined pairs to send and recv
	pair_list->set_hash();

	double dist;
	int rid;
	double n[3];

	for (iwall = 0; iwall < nwalls; iwall++) {
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				itag = tag[i];
				tag2str(str, itag, iwall);
				ipair = pair_list->find_hash(str);
				if (wall_flag[iwall] == BOUND) {
					dim = wall[iwall] / 2;
					side = wall[iwall] % 2;
					for (int j = 0; j < 3; j++) n[j] = 0.0;
					if (side == 0) {
						dist = x[i][dim] - coords0[iwall];
						n[dim] = 1.0;
					}
					else {
						dist = coords0[iwall] - x[i][dim];
						n[dim] = -1.0;
					}
				}
				else if (wall_flag[iwall] == REGION) {
					rid = wall_rid[iwall];
					dist = domain->regions[rid]->find_interaction_distance(n, x[i]);
				}

				// check if the particle is within the interaction range
				if (radius_flag == 0 && fabs(dist) < cut) {
					if (pair_list->nbuilds_total == 1) ipair = pair_list->insert_hash(str);
					else ipair = pair_list->set_pair(str);
				}
				else if (radius_flag == 1 && fabs(dist) < radius[i]) {
					if (pair_list->nbuilds_total == 1 && ipair == -1) ipair = pair_list->insert_hash(str);
					else ipair = pair_list->set_pair(str);
				}
				else {
					if (ipair > -1) {
						pair_list->set_zero(ipair);
					}
				}
			} // if (mask[i] & groupbit)
		} // for (int i = 0; i < nlocal; i++)
	} // for (iwall = 0; iwall < nwalls; iwall++)

	if (pair_list->nbuilds_total > 1 && parallel->nprocs > 1) pair_list->exchange_pair();
}

/* ----------------------------------------------------------------------
   str = itag-fix%name-Wall%id-Side%
------------------------------------------------------------------------- */

void FixWallDEM::tag2str(char *str, int itag, int iwall)
{
	sprintf(str, "%d-fix%s-Wall%d", itag, name, iwall);
}

/* ----------------------------------------------------------------------
   dem: force
------------------------------------------------------------------------- */

void FixWallDEM::pre_force_dem_lsd(int ipair, int iparticle, double delta, double *n, int iwall)
{
	int table, index;

	double drijn, drijnx, drijny, drijnz;                // normal displacement
	double drijtx, drijty, drijtz;                       // tangential displacement
	double vijx, vijy, vijz;                             // relative velocity: vij = vj - vi
	double vijn, vijnx, vijny, vijnz;                    // relative velocity along normal direction
	double vijt, vijt_inv, vijtx, vijty, vijtz;          // relative velocity along tangential direction
	double fijn, fijnx, fijny, fijnz;                    // normal force
	double fijt_a[3], fijt, fijtx, fijty, fijtz;         // tangential force	
	double nx, ny, nz;                                   // unit vector along normal direction
	double tx, ty, tz;                                   // unit vector along tangential direction
	double omegainij[3], omegajnij[3];                   // omega_i cross nij
	double torqueij[3];                                  // torque
	double Li[3], Lj[3];                                 // vector from the center of particle to the contact point 

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double *radius = particle->radius;
	double **omega = particle->omega;
	double **torque = particle->torque;

	double dt = update->dt;

	Region **regions = domain->regions;

	PairList::HashPair **hash_pair = pair_list->hash_pair;
	table = TABLE(ipair);
	index = INDEX(ipair);

	nx = n[0];
	ny = n[1];
	nz = n[2];

	// Overlap
	drijn = delta;
	drijnx = drijn * nx;
	drijny = drijn * ny;
	drijnz = drijn * nz;
	
	// may need to consider wall translational velocity in the future
	double vel_wall[3];
	for (int i = 0; i < 3; i++) vel_wall[i] = 0.0;

	if (wall_flag[iwall] == REGION) {
		int rid = wall_rid[iwall];
		if (regions[rid]->dynamic_flag) {
			for (int j = 0; j < 3; j++) vel_wall[j] = regions[rid]->v_coords[0][j];
		}
	}
	// vij = vi - vj + (wi cross Ri) - (wj cross Rj) (relative velocity)
	vijx = v[iparticle][0] - vel_wall[0];
	vijy = v[iparticle][1] - vel_wall[1];
	vijz = v[iparticle][2] - vel_wall[2];

	if (particle->radius_flag) {
		for (int i = 0; i < 3; i++) Li[i] = (radius[iparticle] - drijn)*(-n[i]);
		Vec_Cross_Prod_3D(omegainij, omega[iparticle], Li);
		vijx += omegainij[0];
		vijy += omegainij[1];
		vijz += omegainij[2];
		if (wall_flag[iwall] == REGION) {
			int rid = wall_rid[iwall];
			if (regions[rid]->rotate_flag) {
				double pro_coord[3], nn[3];
				regions[rid]->find_projection_point(pro_coord, nn, x[iparticle]);
				regions[rid]->rot_axis->find_vector(Lj, pro_coord);
				Vec_Cross_Prod_3D(omegajnij, regions[rid]->omega, Lj);
				// the direction of the center of the region to the contact point is opposite of n
				vijx += -omegajnij[0];
				vijy += -omegajnij[1];
				vijz += -omegajnij[2];
			}
		}
	}

	// |vijn| = vij . nij
	vijn = vijx * nx + vijy * ny + vijz * nz;
    // vijn = |vijn| . nij
	vijnx = vijn*nx;
	vijny = vijn*ny;
	vijnz = vijn*nz;

	if (iparticle == 0 && update->ntimestep > 172700 && iwall == 5) {
		iparticle = iparticle;
	}

	// vijt = vij - (vij . nij) . nij
	vijtx = vijx - vijnx;
	vijty = vijy - vijny;
	vijtz = vijz - vijnz;
	vijt = sqrt(vijtx*vijtx + vijty*vijty + vijtz*vijtz);
	if (fabs(vijt) < EPSILON) vijt_inv = 0.0;
	else vijt_inv = 1.0/vijt;

	// tij = tij/|tij|
	tx = vijtx * vijt_inv;
	ty = vijty * vijt_inv;
	tz = vijtz * vijt_inv;
		
	// calculate drijt
	// first time to contact
	int contact_flag = hash_pair[table][index].ivalues[0];
	double temp;
	if (contact_flag == 0) {
		if (fabs(vijn) < EPSILON) temp = dt;
		else temp = MIN(fabs(drijn/vijn),dt);
		drijtx = vijtx * temp;
		drijty = vijty * temp;
		drijtz = vijtz * temp;
		hash_pair[table][index].ivalues[0] = 1;
	}
	// during the same contact
	else {
		drijtx = hash_pair[table][index].dvalues[0];
		drijty = hash_pair[table][index].dvalues[1];
		drijtz = hash_pair[table][index].dvalues[2];
		
		// update the tang. disp. for the next time step
		drijtx = drijtx + vijtx*dt;
		drijty = drijty + vijty*dt;
		drijtz = drijtz + vijtz*dt;

		temp = drijtx*nx + drijty*ny + drijtz*nz;

		drijtx = drijtx - temp*nx;
		drijty = drijty - temp*ny;
		drijtz = drijtz - temp*nz;
	}

	hash_pair[table][index].dvalues[0] = drijtx;
	hash_pair[table][index].dvalues[1] = drijty;
	hash_pair[table][index].dvalues[2] = drijtz;
	
	// fijn = kn*delta_rijn - Cn*vijn
	fijnx = kn*drijnx - Cn*vijnx;
	fijny = kn*drijny - Cn*vijny;
	fijnz = kn*drijnz - Cn*vijnz;

	fijn = sqrt(fijnx*fijnx + fijny*fijny + fijnz*fijnz);
	
	// tangential force
	// fijt = -kt*delta_rijt - Ct*vijt
	fijtx = -kt * drijtx - Ct * vijtx;
	fijty = -kt * drijty - Ct * vijty;
	fijtz = -kt * drijtz - Ct * vijtz;

	fijt = sqrt(fijtx*fijtx + fijty*fijty + fijtz*fijtz); 
	temp = fabs(mu*fijn);
	// |fijt| > mu*|fnij|
	if (fijt > temp) {
		fijtx = -temp * tx;
		fijty = -temp * ty;
		fijtz = -temp * tz;
	}

	f[iparticle][0] += (fijnx + fijtx);
	f[iparticle][1] += (fijny + fijty);
	f[iparticle][2] += (fijnz + fijtz);

	if (particle->torque_flag) {
		fijt_a[0] = fijtx;
		fijt_a[1] = fijty;
		fijt_a[2] = fijtz;
		Vec_Cross_Prod_3D(torqueij, Li, fijt_a);

		torque[iparticle][0] += torqueij[0];
		torque[iparticle][1] += torqueij[1];
		torque[iparticle][2] += torqueij[2];
	}
}

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
#include "region.h"
#include "region_line.h"
#include "region_polygon.h"
#include "region_point.h"
#include "update.h"
#include "neighbor.h"
#include "modify.h"

#include "particle.h"

#include "pdps_cuda.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "device_functions.h"

using namespace PDPS_NS;
using namespace PhyConst;
using namespace PsMath_NS;

#define EPSILON 1.0e-6

enum{NONE, ROT, TRA_ACC_SIN};

/* ---------------------------------------------------------------------- */

RegionPolygon::RegionPolygon(PDPS *ps, int narg, char **arg) : Region(ps, narg, arg)
{
	//if ((narg - 2) % 4 != 0) error->all(FLERR, "Illegal polygon points");
	nps_initial = narg - 2;

	point_flag = 0;
	line_flag = 0;
	plane_flag = 1;
	volume_flag = 0;

	edges = NULL;

	coords = new double *[nps_initial];
	v_coords = new double*[nps_initial];
	rot_coords = new double *[nps_initial];
	edges = new RegionLine*[nps_initial];
	for (int i = 0; i < nps_initial; i++) {
		coords[i] = new double[3];
		v_coords[i] = new double[3];
		rot_coords[i] = new double[3];
	}

	// parse options
	int iarg = 2;
	nps = 0;
	int pid;
	region_id = domain->nregions;
	while (iarg < narg) {
		if (!strcmp(arg[iarg], "rot")) {
			rotate_flag = 1;
			dynamic_flag = ROT;
			int rid = domain->find_region(arg[iarg+1]);
			if (rid == -1) error->all(FLERR, "Illegal region id");
			rot_axis = domain->regions[rid];
			omega_target = atof(arg[iarg+2]);
			romega = omega_target;
			iarg += 3;
			while (iarg < narg) {
				if (!strcmp(arg[iarg], "start")) {
					start_flag = 1;
					start = atoi(arg[iarg + 1]);
					stable_start = atoi(arg[iarg+2]);
					omega_target = romega;
					romega = 0.0;
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
		}
		else {
			parse_point(arg[iarg], coords[nps]);
			iarg++;
			nps++;
		}
	}

	// create edges
	int narg_line = narg - nps + 2;
	for (int i = 0; i < nps; i++) {
		char **arg_line = new char*[narg_line];
		for (int j = 0; j < narg_line; j++) arg_line[j] = new char[1024];
		sprintf(arg_line[0], "%s-polygon-edges-%d", arg[0], nps);
		sprintf(arg_line[1], "line");
		strcpy(arg_line[2], arg[i%4+2]);
		strcpy(arg_line[3], arg[(i+1)%4+2]);
		for (int j = 4; j < narg_line; j++) {
			strcpy(arg_line[j], arg[j-2+nps]);
		}
		edges[i] = new RegionLine(ps, narg_line, arg_line);
		for (int j = 0; j < narg_line; j++) delete[] arg_line[j];
		delete[] arg_line;
	}


	// initialize translation and rotation velocity
	for (int i = 0; i < nps; i++)
	for (int j = 0; j < 3; j++) {
		v_coords[i][j] = 0.0;
	}

	// calculate normal direction
	calc_normal();
	find_extent_bound();

	// find area
	area = calc_area();

	if (dynamic_flag == ROT) {
		if (start_flag) {
			if (start > stable_start) error->all(FLERR, "Illegal region clylinder rotating start options");
		}
		if (stop_flag) {
			if (stable_end > stop) error->all(FLERR, "Illegal region clylinder rotating stop options");
		}
	}

	char str1[128];
	sprintf(str1, "Its normal direction: %f %f %f\n\n", a, b, c);
	output->print(str1);

	// check if all points are in the same plane ax + by + cz = d
	double temp_d;
	// points (i = 0, 1, 2) actually should already be on the plane
	for (int i = 1; i < nps; i++) {
		temp_d = a*coords[i][0] + b*coords[i][1] + c*coords[i][2];
		if (fabs(temp_d - d) > EPSILON) {
			char str2[128];
			sprintf(str2, "Point %d (%f %f %f) is not on this plane", i+1, coords[i][0], coords[i][1], coords[i][2]);
			error->all(FLERR, str2);
		}
	}

	// find the 2D projection of the 3D plane
	// currently, it is only used for creating particle
	// In the future, this function may be only invoked when creating particle
	// or for some other purpose. 
	plane3D_map2_plane2D(rot_coords);

	// test = fopen("tran_analyze.txt", "w");
}

/* ---------------------------------------------------------------------- */

RegionPolygon::~RegionPolygon()
{
	for (int i = 0; i < nps_initial; i++) {
		delete[] coords[i];
		coords[i] = NULL;
		delete[] rot_coords[i];
		rot_coords[i] = NULL;
	}
	for (int i = 0; i < nps; i++) {
		delete edges[i];
		edges[i] = NULL;
	}
	delete[] coords;
	coords = NULL;
	delete[] rot_coords;
	rot_coords = NULL;
	delete[] v_coords;
	v_coords = NULL;
	delete[] edges;
	edges = NULL;
}

__global__ void gpu_rotate_point_around_axis(
	const int itype, const int rid, const double theta, const int nlocal,
	double *devCoordX, double *devCoordY, double *devCoordZ, int *devType,
	double *devCoord101X, double *devCoord101Y, double *devCoord101Z,
	double *devAxisNormX, double *devAxisNormY, double *devAxisNormZ){

	double x, y, z;
	double cost, sint;

	float a = devCoord101X[rid];
	float b = devCoord101Y[rid];
	float c = devCoord101Z[rid];
	float u = devAxisNormX[rid];
	float v = devAxisNormY[rid];
	float w = devAxisNormZ[rid];

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	for (i = i; i < nlocal; i += blockDim.x * gridDim.x){
		if (devType[i] == itype){
			x = devCoordX[i];
			y = devCoordY[i];
			z = devCoordZ[i];

			cost = cos(theta);
			sint = sin(theta);

			devCoordX[i] = (a*(v*v + w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1 - cost) +
				x*cost + (-c*v + b*w - w*y + v*z)*sint;
			devCoordY[i] = (b*(u*u + w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1 - cost) +
				y*cost + (c*u - a*w + w*x - u*z)*sint;
			devCoordZ[i] = (c*(u*u + v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1 - cost) +
				z*cost + (-b*u + a*v - v*x + u*y)*sint;
		}

	}

}

__device__ void gpu_calc_normal(
	int rid,
	double *devCoord1X, double *devCoord1Y, double *devCoord1Z,
	double *devCoord2X, double *devCoord2Y, double *devCoord2Z,
	double *devCoord3X, double *devCoord3Y, double *devCoord3Z,
	double *devCoord4X, double *devCoord4Y, double *devCoord4Z,
	double *devA, double *devB, double *devC, double *devD){

	double n12[3], n13[3], normal[3];
	n12[0] = devCoord3X[rid] - devCoord2X[rid];
	n12[1] = devCoord3Y[rid] - devCoord2Y[rid];
	n12[2] = devCoord3Z[rid] - devCoord2Z[rid];
	n13[0] = devCoord4X[rid] - devCoord2X[rid];
	n13[1] = devCoord4Y[rid] - devCoord2Y[rid];
	n13[2] = devCoord4Z[rid] - devCoord2Z[rid];

	normal[0] = n12[1] * n13[2] - n12[2] * n13[1];
	normal[1] = n12[2] * n13[0] - n12[0] * n13[2];
	normal[2] = n12[0] * n13[1] - n12[1] * n13[0];

	double sq = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

	normal[0] /= sq;
	normal[1] /= sq;
	normal[2] /= sq;

	// plane: ax + by + cz = d
	devA[rid] = normal[0];
	devB[rid] = normal[1];
	devC[rid] = normal[2];
	devD[rid] = devA[rid] * devCoord1X[rid] + devB[rid] * devCoord1Y[rid] + devC[rid] * devCoord1Z[rid];
}

__global__ void gpuDynamicCheck(
	double *devCoordX, double *devCoordY, double *devCoordZ,
	double *devVeloX, double *devVeloY, double *devVeloZ,
	int *devWallContactFlag, double *devWallDrijtx, double *devWallDrijty, double *devWallDrijtz,
	int *devWall, double *devCoords0, int *devWall_rid, int *devWall_flag,
	double *devMass, double *devRadius, int *devType, int *devMask, const int nlocal, const double dt,
	const int nwalls,
	double *devForceX, double *devForceY, double *devForceZ,
	int *devStyle, double *devRadiusCylinder, double *devHeight,
	double *devCoord1X, double *devCoord1Y, double *devCoord1Z,
	double *devCoord2X, double *devCoord2Y, double *devCoord2Z,
	double *devCoord3X, double *devCoord3Y, double *devCoord3Z,
	double *devCoord4X, double *devCoord4Y, double *devCoord4Z,
	double *devA, double *devB, double *devC, double *devD,
	const int ntimestep,
	int *devRotateFlag, int *devStartFlag, int *devEndFlag,
	int *devStart, int *devStableStart, int *devEnd, int *devStableEnd,
	double *devOmegaTarget,
	double *devCoord101X, double *devCoord101Y, double *devCoord101Z,
	double *devCoord102X, double *devCoord102Y, double *devCoord102Z,
	double *devAxisNormX, double *devAxisNormY, double *devAxisNormZ){
	int iwall = blockIdx.x * blockDim.x + threadIdx.x;

	//double cos_beta, sin_beta;
	//double dist, dt;
	double theta;
	double romega, domega, omega[3], normal[3];
	int ALI1, ALI2, ALI3;
	double ALI11, ALI12, ALI13;
	if (iwall < nwalls) {
		ALI1 = iwall;
		int rid = devWall_rid[iwall];
		ALI2 = rid;
		normal[0] = devAxisNormX[rid];
		normal[1] = devAxisNormY[rid];
		normal[2] = devAxisNormZ[rid];
		double omega_target = devOmegaTarget[rid];
		double start = devStart[rid];
		double stable_start = devStableStart[rid];
		double stop = devEnd[rid];
		double stable_end = devStableEnd[rid];

		if (devWall_flag[iwall] == 1 && devStyle[rid] == 2 &&
			devRotateFlag[rid] == 1) {
			//		if (dynamic_flag == ROT) {
			// update rotating velocity
			if (devStartFlag[rid] == 0) romega = omega_target;
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
			if (devEndFlag[rid]) {
				if (ntimestep >= stop) {
					romega = 0.0;
					//dynamic_flag = NONE;
				}
				else if (ntimestep >= stable_end && ntimestep < stop) {
					domega = omega_target / (stop - stable_end);
					romega = omega_target - domega * (ntimestep - stable_end);
				}
			}

			for (int i = 0; i < 3; i++) omega[i] = romega*normal[i];

			// rotate theta angle
			theta = romega * dt;
			//gpu_rotate_point_around_axis(1, 0, rid, theta, devCoordX, devCoordY, devCoordZ,
			//	devCoord1X, devCoord1Y, devCoord1Z, devCoord101X, devCoord101Y, devCoord101Z,
			//	devCoord102X, devCoord102Y, devCoord102Z, devAxisNormX, devAxisNormY, devAxisNormZ);
			//gpu_rotate_point_around_axis(1, 0, rid, theta, devCoordX, devCoordY, devCoordZ,
			//	devCoord2X, devCoord2Y, devCoord2Z, devCoord101X, devCoord101Y, devCoord101Z,
			//	devCoord102X, devCoord102Y, devCoord102Z, devAxisNormX, devAxisNormY, devAxisNormZ);
			//gpu_rotate_point_around_axis(1, 0, rid, theta, devCoordX, devCoordY, devCoordZ,
			//	devCoord3X, devCoord3Y, devCoord3Z, devCoord101X, devCoord101Y, devCoord101Z,
			//	devCoord102X, devCoord102Y, devCoord102Z, devAxisNormX, devAxisNormY, devAxisNormZ);
			//gpu_rotate_point_around_axis(1, 0, rid, theta, devCoordX, devCoordY, devCoordZ,
			//	devCoord4X, devCoord4Y, devCoord4Z, devCoord101X, devCoord101Y, devCoord101Z,
			//	devCoord102X, devCoord102Y, devCoord102Z, devAxisNormX, devAxisNormY, devAxisNormZ);

			// This part is used to visualize the blades.
			// I did not delete it is because we may need it to visualize sth for debug use

			for (int i = 0; i < nlocal; i++) {
				if (devType[i] == 3) {
			/*		gpu_rotate_point_around_axis(2, i, rid, theta, devCoordX, devCoordY, devCoordZ,
						devCoord1X, devCoord1Y, devCoord1Z, devCoord101X, devCoord101Y, devCoord101Z,
						devCoord102X, devCoord102Y, devCoord102Z, devAxisNormX, devAxisNormY, devAxisNormZ);*/
				}
			}
			gpu_calc_normal(rid,
				devCoord1X, devCoord1Y, devCoord1Z,
				devCoord2X, devCoord2Y, devCoord2Z,
				devCoord3X, devCoord3Y, devCoord3Z,
				devCoord4X, devCoord4Y, devCoord4Z,
				devA, devB, devC, devD);
			int ALI99 = 100;
		}
	}
	int ALI100 = ALI1 + ALI2 + ALI3;
	double ALI110 = ALI11 + ALI12 + ALI13;
	int ALI200 = ALI100 / ALI110;
	int ALI300 = ALI200 + ALI100 + ALI110;
	if (ALI1>0 && ALI200>50) ALI200 = ALI300;
}

/* ----------------------------------------------------------------------
   Dynamic check: 
   update rotation speed
   update rotation angle
------------------------------------------------------------------------- */

void RegionPolygon::dynamic_check()
{
	double cos_beta, sin_beta;
	double dist, dt;
	double theta;

	if (dynamic_flag == ROT) {
		// update rotating velocity
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


		// This part is used to visualize the blades.
		// I did not delete it is because we may need it to visualize sth for debug use

		//for (int i = 0; i < particle->nlocal; i++) {
		//	if (particle->type[i] == 3 && strcmp(name, "blade1") == 0) {
		//		rot_axis->rotate_point_around_axis(particle->x[i], theta);
		//	}
		//}
		cudaError_t error_t;
		gpu_rotate_point_around_axis << < int(particle->nlocal + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >
			(3, region_id, theta, particle->nlocal,
			particle->devCoordX, particle->devCoordY, particle->devCoordZ, particle->devType,
			domain->devCoord101X, domain->devCoord101Y, domain->devCoord101Z,
			domain->devAxisNormX, domain->devAxisNormY, domain->devAxisNormZ);
		error_t = cudaMemcpy(neighbor->hostCoordX, particle->devCoordX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
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
			cur_time += dt;
		}
	}

	// update the edges
	for (int i = 0; i < nps; i++) {
		edges[i]->dynamic_check();
	}

	calc_normal();
	find_extent_bound();


	//	update relavent parameters in GPU
	domain->hostCoord1X[region_id] = domain->regions[region_id]->coords[0][0];
	domain->hostCoord1Y[region_id] = domain->regions[region_id]->coords[0][1];
	domain->hostCoord1Z[region_id] = domain->regions[region_id]->coords[0][2];
	domain->hostCoord2X[region_id] = domain->regions[region_id]->coords[1][0];
	domain->hostCoord2Y[region_id] = domain->regions[region_id]->coords[1][1];
	domain->hostCoord2Z[region_id] = domain->regions[region_id]->coords[1][2];
	domain->hostCoord3X[region_id] = domain->regions[region_id]->coords[2][0];
	domain->hostCoord3Y[region_id] = domain->regions[region_id]->coords[2][1];
	domain->hostCoord3Z[region_id] = domain->regions[region_id]->coords[2][2];
	domain->hostCoord4X[region_id] = domain->regions[region_id]->coords[3][0];
	domain->hostCoord4Y[region_id] = domain->regions[region_id]->coords[3][1];
	domain->hostCoord4Z[region_id] = domain->regions[region_id]->coords[3][2];
	domain->hostVelo0X[region_id] = domain->regions[region_id]->v_coords[0][0];
	domain->hostVelo0Y[region_id] = domain->regions[region_id]->v_coords[0][1];
	domain->hostVelo0Z[region_id] = domain->regions[region_id]->v_coords[0][2];
	domain->hostA[region_id] = domain->regions[region_id]->a;
	domain->hostB[region_id] = domain->regions[region_id]->b;
	domain->hostC[region_id] = domain->regions[region_id]->c;
	domain->hostD[region_id] = domain->regions[region_id]->d;

	cudaMemcpy(domain->devCoord1X, domain->hostCoord1X, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord1Y, domain->hostCoord1Y, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord1Z, domain->hostCoord1Z, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord2X, domain->hostCoord2X, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord2Y, domain->hostCoord2Y, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord2Z, domain->hostCoord2Z, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord3X, domain->hostCoord3X, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord3Y, domain->hostCoord3Y, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord3Z, domain->hostCoord3Z, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord4X, domain->hostCoord4X, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord4Y, domain->hostCoord4Y, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devCoord4Z, domain->hostCoord4Z, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devVelo0X, domain->hostVelo0X, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devVelo0Y, domain->hostVelo0Y, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devVelo0Z, domain->hostVelo0Z, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devA, domain->hostA, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devB, domain->hostB, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devC, domain->hostC, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(domain->devD, domain->hostD, domain->nregions * sizeof(double), cudaMemcpyHostToDevice);
	//gpuDynamicCheck << < int(particle->nlocal + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >
	//	(particle->devCoordX, particle->devCoordY, particle->devCoordZ,
	//	particle->devVeloX, particle->devVeloY, particle->devVeloZ,
	//	neighbor->devWallContactFlag, neighbor->devWallDrijtx, neighbor->devWallDrijty, neighbor->devWallDrijtz,
	//	neighbor->devWall, neighbor->devCoords0, neighbor->devWall_rid, neighbor->devWall_flag,
	//	particle->devMass, particle->devRadius, particle->devType,
	//	particle->devMask, particle->nlocal, update->dt,
	//	modify->gpunwalls,
	//	particle->devForceX, particle->devForceY, particle->devForceZ,
	//	domain->devStyle, domain->devRadiusCylinder, domain->devHeight,
	//	domain->devCoord1X, domain->devCoord1Y, domain->devCoord1Z,
	//	domain->devCoord2X, domain->devCoord2Y, domain->devCoord2Z,
	//	domain->devCoord3X, domain->devCoord3Y, domain->devCoord3Z,
	//	domain->devCoord4X, domain->devCoord4Y, domain->devCoord4Z,
	//	domain->devA, domain->devB, domain->devC, domain->devD,
	//	update->ntimestep,
	//	domain->devRotateFlag, domain->devStartFlag, domain->devEndFlag,
	//	domain->devStart, domain->devStableStart, domain->devEnd, domain->devStableEnd,
	//	domain->devOmegaTarget,
	//	domain->devCoord101X, domain->devCoord101Y, domain->devCoord101Z,
	//	domain->devCoord102X, domain->devCoord102Y, domain->devCoord102Z,
	//	domain->devAxisNormX, domain->devAxisNormY, domain->devAxisNormZ);
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegionPolygon::projection_inside(double *x0)
{
	int flag, side_flag;

	double pro_coord[3], n[3];

	// side_flag = 1: in the space where the normal direction points;
	//             0: in the other side
	find_projection_point(pro_coord, n, x0);

	flag = inside(pro_coord);

	return flag;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegionPolygon::projection_inside(double *x0, double radius)
{
	int flag, side_flag;

	flag = 0;
	// check the edges first
	double dist;
	for (int i = 0; i < nps; i++) {
		dist = edges[i]->find_distance(x0);
		if (fabs(dist) < radius) {
			return 1;
		}
	}

	double pro_coord[3], n[3];

	// side_flag = 1: in the space where the normal direction points;
	//             0: in the other side
	find_projection_point(pro_coord, n, x0);

	flag = inside(pro_coord);

	return flag;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegionPolygon::inside(double *x0)
{
	int i, j;
	int flag;

	double p1_norm, p2_norm, p1p2;
	double p1[3], p2[3];
	double c[3];
	double cos_theta, theta;
	theta = 0.0;

	for (i = 0; i < nps; i++) {
		for (int k = 0; k < 3; k++) p1[k] = coords[i][k] - x0[k];

		j = (i+1) % nps;
		for (int k = 0; k < 3; k++) p2[k] = coords[j][k] - x0[k];

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

void RegionPolygon::calc_normal()
{
	double n12[3], n13[3];
	for (int i = 0; i < 3; i++) {
		n12[i] = coords[2][i] - coords[1][i];
		n13[i] = coords[3][i] - coords[1][i];
	}
	Vec_Cross_Prod_3D(normal, n12, n13);

	double sq = Vec_Norm2(normal);

	if (sq < EPSILON) error->all(FLERR, "Points coordinates are not correct, they meight be too close to each other");

	normal[0] /= sq;
	normal[1] /= sq;
	normal[2] /= sq;

	// plane: ax + by + cz = d
	a = normal[0];
	b = normal[1];
	c = normal[2];
	d = a*coords[0][0] + b*coords[0][1] + c*coords[0][2];
}

/* ---------------------------------------------------------------------- */

double RegionPolygon::calc_area()
{
	double s;
	double pij[3], temp[3];

	s = 0.0;

	for (int i = 0; i < 3; i++) pij[i] = 0.0;
	int j = 0;
	for (int i = 0; i < nps; i++) {
		j = (i+1) % nps;
		Vec_Cross_Prod_3D(temp, coords[j], coords[i]);
		for (int k = 0; k < 3; k++) pij[k] += temp[k];
	}

	s = Vec_Dot_Prod_3D(normal, pij);
	s = fabs(s)/2;

	return s;
}

/* ---------------------------------------------------------------------- */

void RegionPolygon::find_extent_bound()
{
	for (int i = 0; i < 3; i++) {
		extent_lo[i] = coords[0][i];
		extent_hi[i] = coords[0][i];
	}

	// outbound of the region
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

/* ---------------------------------------------------------------------- */

double RegionPolygon::find_distance(double *x0)
{
	double dist;

	dist = a*x0[0] + b*x0[1] + c*x0[2] - d; 

	return dist;
}

/* ---------------------------------------------------------------------- */

double RegionPolygon::find_interaction_distance(double *n, double *x0)
{
	int flag;
	double dist, pro_coord[3];

	dist = find_projection_point(pro_coord, n, x0);
	flag = inside(pro_coord);

	if (flag == 0) {
		double n1[3];
		// check the edges first
		dist = edges[0]->find_interaction_distance(n, x0);
		double temp;
		for (int i = 1; i < nps; i++) {
			temp = edges[i]->find_interaction_distance(n1, x0);
			if (fabs(dist) > fabs(temp)) {
				dist = temp;
				for (int j = 0; j < 3; j++) n[j] = n1[j];
			}
		}
	}
	
	return dist;
}

/* ---------------------------------------------------------------------- */

void RegionPolygon::find_vector_to_rot_axis(double *vec, double *x0)
{
	rot_axis->find_vector(vec, x0);
}

/* ---------------------------------------------------------------------- */

double RegionPolygon::find_projection_point(double *p, double *n, double *x0)
{
	double dist[3], dist_norm;

	dist_norm = find_distance(x0);

	for (int i = 0; i < 3; i++) {
		n[i] = dist_norm / fabs(dist_norm) * normal[i];
		dist[i] = dist_norm*normal[i];
		p[i] = x0[i] - dist[i];
	}
	return dist_norm;
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

void RegionPolygon::plane3D_map2_plane2D(double **p_coords)
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

	temp1 = 0.0;
	for (i = 0; i < 3; i++) {
		temp1 += (coords[1][i] - coords[0][i]) * (coords[1][i] - coords[0][i]);
	}
	temp1 = sqrt(temp1);
	inv_temp1 = 1.0 / temp1;

	for (i = 0; i < 3; i++) temp2[i] = (coords[1][i] - coords[0][i]) * inv_temp1;

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

/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include <iostream>

#include "domain.h"
#include "error.h"
#include "force.h"
#include "integrate.h"
#include "neighbor.h"
#include "region.h"
#include "update.h"
#include "particle.h"
#include "modify.h"
#include "v_verlet.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

Update::Update(PDPS *ps) : Pointers(ps)
{
	ntimestep = 0;

	integrate_style = NULL;
	integrate = NULL;
	//integrate = new Integrate(ps);

	units_flag = 0;
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
	delete integrate_style;
	delete integrate;
	integrate = NULL;
	integrate_style = NULL;
}

/* ---------------------------------------------------------------------- */

void Update::init()
{
	integrate->init();
}

/* ----------------------------------------------------------------------
   Dynamic check: 
   update rotation speed
   update rotation angle
------------------------------------------------------------------------- */
__device__ void gpu_rotate_point_around_axis(
	int flag, int i, int rid, double theta,
    double *devCoordX, double *devCoordY, double *devCoordZ,
	double *devCoord1X, double *devCoord1Y, double *devCoord1Z,
	double *devCoord101X, double *devCoord101Y, double *devCoord101Z, 
	double *devCoord102X, double *devCoord102Y, double *devCoord102Z, 
	double *devAxisNormX, double *devAxisNormY, double *devAxisNormZ){
	
	double a, b, c;
	double u, v, w;
	double x, y, z;
	double cost, sint;
    
	a = devCoord101X[rid];
	b = devCoord101Y[rid];
	c = devCoord101Z[rid];
	u = devAxisNormX[rid];
	v = devAxisNormY[rid];
	w = devAxisNormZ[rid];
    
	if (flag == 1) {
		x = devCoord1X[rid];
		y = devCoord1Y[rid];
		z = devCoord1Z[rid];
	} else {
		x = devCoordX[i];
		y = devCoordY[i];
		z = devCoordZ[i];
	}
    
	cost = cos(theta);
	sint = sin(theta);
	
	if (flag == 1) {
		devCoord1X[rid] = (a*(v*v + w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1 - cost) +
							x*cost + (-c*v + b*w - w*y + v*z)*sint;
		devCoord1Y[rid] = (b*(u*u + w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1 - cost) +
							y*cost + (c*u - a*w + w*x - u*z)*sint;
		devCoord1Z[rid] = (c*(u*u + v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1 - cost) +
							z*cost + (-b*u + a*v - v*x + u*y)*sint;
	} else {
		devCoordX[i] = (a*(v*v + w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1 - cost) +
						x*cost + (-c*v + b*w - w*y + v*z)*sint;
		devCoordY[i] = (b*(u*u + w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1 - cost) +
						y*cost + (c*u - a*w + w*x - u*z)*sint;
		devCoordZ[i] = (c*(u*u + v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1 - cost) +
						z*cost + (-b*u + a*v - v*x + u*y)*sint;
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
	
	normal[0] = n12[1]*n13[2] - n12[2]*n13[1];
	normal[1] = n12[2]*n13[0] - n12[0]*n13[2];
	normal[2] = n12[0]*n13[1] - n12[1]*n13[0];

	double sq = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);

	normal[0] /= sq;
	normal[1] /= sq;
	normal[2] /= sq;

	// plane: ax + by + cz = d
	devA[rid] = normal[0];
	devB[rid] = normal[1];
	devC[rid] = normal[2];
	devD[rid] = devA[rid]*devCoord1X[rid] + devB[rid]*devCoord1Y[rid] + devC[rid]*devCoord1Z[rid];
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
		ALI1=iwall;
		int rid = devWall_rid[iwall];
		ALI2=rid;
		normal[0] = devAxisNormX[rid];
		normal[1] = devAxisNormY[rid];
		normal[2] = devAxisNormZ[rid];
		double omega_target = devOmegaTarget[rid];
		double start = devStart[rid];
		double stable_start = devStableStart[rid];
		double stop = devEnd[rid];
		double stable_end = devStableEnd[rid];
		
		if (devWall_flag[iwall] == 1 && devStyle[rid]==2 &&
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
			gpu_rotate_point_around_axis(1, 0, rid, theta, devCoordX, devCoordY, devCoordZ,
										devCoord1X, devCoord1Y, devCoord1Z, devCoord101X, devCoord101Y, devCoord101Z,
										devCoord102X, devCoord102Y, devCoord102Z, devAxisNormX, devAxisNormY, devAxisNormZ);
			gpu_rotate_point_around_axis(1, 0, rid, theta, devCoordX, devCoordY, devCoordZ,
										devCoord2X, devCoord2Y, devCoord2Z, devCoord101X, devCoord101Y, devCoord101Z,
										devCoord102X, devCoord102Y, devCoord102Z, devAxisNormX, devAxisNormY, devAxisNormZ);
			gpu_rotate_point_around_axis(1, 0, rid, theta, devCoordX, devCoordY, devCoordZ,
										devCoord3X, devCoord3Y, devCoord3Z, devCoord101X, devCoord101Y, devCoord101Z,
										devCoord102X, devCoord102Y, devCoord102Z, devAxisNormX, devAxisNormY, devAxisNormZ);
			gpu_rotate_point_around_axis(1, 0, rid, theta, devCoordX, devCoordY, devCoordZ,
										devCoord4X, devCoord4Y, devCoord4Z, devCoord101X, devCoord101Y, devCoord101Z,
										devCoord102X, devCoord102Y, devCoord102Z, devAxisNormX, devAxisNormY, devAxisNormZ);
	
			// This part is used to visualize the blades.
			// I did not delete it is because we may need it to visualize sth for debug use
	
			for (int i = 0; i < nlocal; i++) {
				if (devType[i] == 3) {
					gpu_rotate_point_around_axis(2, i, rid, theta, devCoordX, devCoordY, devCoordZ,
										devCoord1X, devCoord1Y, devCoord1Z, devCoord101X, devCoord101Y, devCoord101Z,
										devCoord102X, devCoord102Y, devCoord102Z, devAxisNormX, devAxisNormY, devAxisNormZ);
				}
			}
//		}
//		else if (dynamic_flag == TRA_ACC_SIN) {
//			double acc[3];
//			double dt;
//	
//			if (stop_flag) {
//				if (update->ntimestep >= stop) {
//					dynamic_flag = NONE;
//				}
//			}
//	
//			if (update->ntimestep >= start && ((stop_flag == 1 && update->ntimestep <= stop) || stop_flag == 0)) {
//				dt = update->dt;
//				for (int i = 0; i < 3; i++) {
//					acc[i] = tra_A[i] * sin(tra_wt[i] * cur_time + tra_phi[i]);
//				}
//				// translate the coordinate of points
//				for (int i = 0; i < nps; i++)
//				for (int j = 0; j < 3; j++) {
//					v_coords[i][j] += acc[j] * dt;
//					coords[i][j] += v_coords[i][j] * dt;
//				}
//				cur_time += dt;
//			}
//		}
		
		// update the edges
		//for (int i = 0; i < nps; i++) {
			//edges[i]->dynamic_check();
		//}
	
		//calc_normal();
		//find_extent_bound();
			gpu_calc_normal(rid,
							devCoord1X, devCoord1Y, devCoord1Z,
							devCoord2X, devCoord2Y, devCoord2Z,
							devCoord3X, devCoord3Y, devCoord3Z,
							devCoord4X, devCoord4Y, devCoord4Z,
							devA, devB, devC, devD);
			int ALI99 = 100;
		}
	}
	int ALI100 = ALI1+ALI2+ALI3;
	double ALI110 = ALI11+ALI12+ALI13;
	int ALI200 = ALI100/ALI110;
	int ALI300 = ALI200+ALI100+ALI110;
	if (ALI1>0 && ALI200>50) ALI200 = ALI300;
}

/* ---------------------------------------------------------------------- */

void Update::dynamic_check()
{
//	for (int i = 0; i < domain->nregions; i++) {
//		if (domain->regions[i]->dynamic_flag) {
//			domain->regions[i]->dynamic_check();
//		}
//	}
//std::cout<<"Here A .................."<<std::endl;
	gpuDynamicCheck << < int(particle->nlocal + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >
	 (particle->devCoordX, particle->devCoordY, particle->devCoordZ,
	  particle->devVeloX, particle->devVeloY, particle->devVeloZ,
	  neighbor->devWallContactFlag, neighbor->devWallDrijtx, neighbor->devWallDrijty, neighbor->devWallDrijtz,
	  neighbor->devWall, neighbor->devCoords0, neighbor->devWall_rid, neighbor->devWall_flag,
	  particle->devMass, particle->devRadius, particle->devType,
	  particle->devMask, particle->nlocal, update->dt, 
	  modify->gpunwalls,
	  particle->devForceX, particle->devForceY, particle->devForceZ,
	  domain->devStyle, domain->devRadiusCylinder, domain->devHeight,
	  domain->devCoord1X, domain->devCoord1Y, domain->devCoord1Z,
	  domain->devCoord2X, domain->devCoord2Y, domain->devCoord2Z,
	  domain->devCoord3X, domain->devCoord3Y, domain->devCoord3Z,
	  domain->devCoord4X, domain->devCoord4Y, domain->devCoord4Z,
	  domain->devA, domain->devB, domain->devC, domain->devD,
	  ntimestep,
	  domain->devRotateFlag, domain->devStartFlag, domain->devEndFlag,
	  domain->devStart, domain->devStableStart, domain->devEnd, domain->devStableEnd,
	  domain->devOmegaTarget,
	  domain->devCoord101X, domain->devCoord101Y, domain->devCoord101Z, 
	  domain->devCoord102X, domain->devCoord102Y, domain->devCoord102Z, 
	  domain->devAxisNormX, domain->devAxisNormY, domain->devAxisNormZ);
//std::cout<<"Here B .................."<<std::endl;
}

/* ----------------------------------------------------------------------
						  Set units
------------------------------------------------------------------------- */

void Update::set_units(const char *style)
{
	// physical constants from:
	// http://physics.nist.gov/cuu/Constants/Table/allascii.txt
	// using thermochemical calorie = 4.184 J

	if (strcmp(style,"lj") == 0) {
		force->boltz = 1.0;
		//force->boltz = 8.300365358821798e-06;
		force->hplanck = 0.18292026;  // using LJ parameters for argon
		force->mvv2e = 1.0;
		force->ftm2v = 1.0;
		force->mv2d = 1.0;
		force->nktv2p = 1.0;
		force->qqr2e = 1.0;
		force->qe2f = 1.0;
		force->vxmu2f = 1.0;
		force->xxt2kmu = 1.0;
		force->e_mass = 0.0;    // not yet set
		force->hhmrr2e = 0.0;
		force->mvh2r = 0.0;
		force->angstrom = 1.0;
		force->femtosecond = 1.0;
		force->qelectron = 1.0;

		dt = 0.005;
		neighbor->rskin = 0.3;
		units_flag = 1;
	}

	if (units_flag == 0) {
		char str[128];
		sprintf(str, "units style %s is invalid", style);
		error->all(FLERR,str);
	}
}

/* ----------------------------------------------------------------------
						 create integrate class
------------------------------------------------------------------------- */

void Update::create_integrate(int narg, char **arg)
{
	if (narg < 1) error->all(FLERR,"Illegal run_style command");

	delete [] integrate_style;
	integrate_style = NULL;
	delete integrate;
	integrate = NULL;
	
	new_integrate(arg[0], narg, arg);
}

/* ----------------------------------------------------------------------
						create new integrate class
------------------------------------------------------------------------- */

void Update::new_integrate(char *style, int narg, char **arg)
{
	int success = 0;
	
#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) \
    if (strcmp(style,#key) == 0) integrate = new Class(ps,narg,arg);
#include "style_integrate.h"
#undef IntegrateStyle
#undef INTEGRATE_CLASS
	
}

/* ---------------------------------------------------------------------- */

void Update::reset_timestep(int narg, char **arg)
{
	if (narg != 1) error->all(FLERR,"Illegal reset_timestep command");
	int newstep = atoi(arg[0]);
	if (newstep < 0) error->all(FLERR,"Timestep must be >= 0");
	// if new step is too big
	ntimestep = newstep;
}

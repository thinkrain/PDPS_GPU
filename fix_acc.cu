/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "domain.h"
#include "error.h"
#include "fix_acc.h"
#include "particle.h"
#include "region.h"
#include "timer.h"
#include "update.h"
#include "neighbor.h"

#include "pdps_cuda.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
using namespace PDPS_NS;
using namespace FixConst;


__global__ void gpuacc(double *devCoordX, double *devCoordY, double *devCoordZ, int *devMask, const int groupbit,
						const int nlocal, double *devForceX, double *devForceY, double *devForceZ,
						const double xacc, const double yacc, const double zacc, double *devMass, int *devType){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	__shared__ double acc[3];
	__shared__ double mass[10];
	double fx, fy, fz;
	acc[0] = xacc;
	acc[1] = yacc;
	acc[2] = zacc;
	for (int tid = 0; tid < 10; tid++)
		mass[tid] = devMass[tid];
	__syncthreads();
	for (i = i; i < nlocal; i += blockDim.x * gridDim.x){
		if (devMask[i] & groupbit){
			int itype = devType[i];
			fx = mass[itype] * acc[0];
			fy = mass[itype] * acc[1];
			fz = mass[itype] * acc[2];
			devForceX[i] += fx;  
			devForceY[i] += fy;
			devForceZ[i] += fz;
		}

	}


}
/* ---------------------------------------------------------------------- */

FixAcc::FixAcc(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 6) {
		error->all(FLERR,"Illegal fix acc command");
	}
	
	xacc = yacc = zacc = 0.0;

	region_flag = 0;
	atomic_flag = 0;
	sphere_flag = 0;
	int iarg;
	iarg = 3;
	while (iarg < narg) {
		if (!strcmp(arg[iarg],"gravity")) {
			if (!strcmp(arg[iarg+1],"x")) xacc = atof(arg[iarg+2]);
			if (!strcmp(arg[iarg+1],"y")) yacc = atof(arg[iarg+2]);
			if (!strcmp(arg[iarg+1],"z")) zacc = atof(arg[iarg+2]);
			if (narg > iarg + 3){
				if (!strcmp(arg[iarg + 3], "atomic")) atomic_flag = 1;
				if (!strcmp(arg[iarg + 3], "sphere")) sphere_flag = 1;
			}
			iarg += 4;
		}
		else if (!strcmp(arg[iarg], "region")) {
			rid = domain->find_region(arg[iarg+1]);
			if (rid == -1) error->all(FLERR, "Cannot find the region id");
			region_flag = 1;
			iarg += 2;
		}
		else error->all(FLERR, "Illegal command option");
	}

}

/* ---------------------------------------------------------------------- */

FixAcc::~FixAcc()
{

}

/* ---------------------------------------------------------------------- */

int FixAcc::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixAcc::init()
{

}

/* ---------------------------------------------------------------------- */

void FixAcc::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixAcc::post_force()
{
	double **x = particle->x;
	double **f = particle->f;
	double *mass = particle->mass;
	double *rmass = particle->rmass;
	int rmass_flag = particle->rmass_flag;
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	double massone;

	int inside_flag;

	//for (int i = 0; i < nlocal; i++) {
	//	if (mask[i] & groupbit) {
	//		if (region_flag == 1) {
	//			inside_flag = domain->regions[rid]->inside(x[i]);
	//			if (inside_flag == 0) continue;
	//		}
	//		if (sphere_flag) 
	//			massone = rmass[i];
	//		else massone = mass[type[i]];
	//		f[i][0] += massone*xacc;
	//		f[i][1] += massone*yacc;
	//		f[i][2] += massone*zacc;
	//	}
	//}
	//cudaError_t error_t;
	//error_t = cudaMemcpy(neighbor->hostForceX, particle->devForceX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(neighbor->hostForceY, particle->devForceY, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(neighbor->hostForceZ, particle->devForceZ, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	gpuacc << < int(nlocal + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >> >(particle->devCoordX, particle->devCoordY, particle->devCoordZ, particle->devMask, groupbit,
		nlocal, particle->devForceX, particle->devForceY, particle->devForceZ,
		xacc, yacc, zacc, particle->devMass, particle->devType);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	time = time;
	//error_t = cudaMemcpy(neighbor->hostForceX, particle->devForceX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(neighbor->hostForceY, particle->devForceY, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(neighbor->hostForceZ, particle->devForceZ, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//particle->TransferG2C();
}

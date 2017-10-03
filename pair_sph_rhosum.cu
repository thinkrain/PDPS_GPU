/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "parallel.h"
#include "error.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair_sph_rhosum.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;
#include "pdps_cuda.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "device_functions.h"

#define DELTA 1
#define EPSILON 1.0e-10
#define PI 3.1416

__global__ void gpuComputerho(double *devCoordX, double *devCoordY, double *devCoordZ, int *devPairtable, int *devNumneigh,
	double *devRho, double *devMass, int *devType, int *devMask, const double h,
	const int nlocal, const double a3D, int *devSetflag, double *devCutsq){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	double wf, xtemp, ytemp, ztemp, rsq, delx, dely, delz, q;
	int j, jj, itype, jtype, jnum;
	__shared__ double mass[TYPEMAX];
	__shared__ int setflag[TYPEMAX * TYPEMAX];
	__shared__ double cutsq[TYPEMAX * TYPEMAX];
	for (int tid = 0; tid < TYPEMAX; tid++){
		mass[tid] = devMass[tid];
		for (j = 0; j < TYPEMAX; j++){
			setflag[tid * TYPEMAX + j] = devSetflag[tid * TYPEMAX + j];
			cutsq[tid * TYPEMAX + j] = devCutsq[tid * TYPEMAX + j];
		}

	}

	//for (i = i; i < nlocal; i += blockDim.x * gridDim.x){
	//	itype = devType[i];
	//	wf = a3D;
	//	if (setflag[itype * TYPEMAX + itype]){
	//		devRho[i] = mass[itype] * wf;
	//	}
	//		
	//}
	__syncthreads();

	for (i = blockIdx.x * blockDim.x + threadIdx.x; i < nlocal; i += blockDim.x * gridDim.x){

		itype = devType[i];
		if (setflag[itype * TYPEMAX + itype] == 1){
			wf = a3D;
			double irho = mass[itype] * wf;
			irho = mass[itype] * wf;
			xtemp = devCoordX[i];
			ytemp = devCoordY[i];
			ztemp = devCoordZ[i];
			jnum = devNumneigh[i];
			for (jj = 0; jj < jnum; jj++){
				j = devPairtable[i * NEIGHMAX + jj];
				jtype = devType[j];
				if (setflag[itype * TYPEMAX + jtype] == 1){
					delx = xtemp - devCoordX[j];
					dely = ytemp - devCoordY[j];
					delz = ztemp - devCoordZ[j];
					rsq = delx * delx + dely * dely + delz * delz;
					if (rsq < cutsq[itype * TYPEMAX + jtype]){
						q = sqrt(rsq) / h;
						//	Cubic Spline
						if (q < 1)
							wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
						else
							wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
						wf = wf * a3D;
						//  detect solid particle for local average

						irho += mass[jtype] * wf;

					}		//  rsq < cutsq[itype][jtype]

				}	// setflag[itype * 10 + jtype]
			}	// j < jnum
			/*if (irho < 950)
				irho = 1000;*/
			devRho[i] = irho;
		}	//	setflag[itype * TYPEMAX + itype]

	}	// i < nlocal

}
/* ---------------------------------------------------------------------- */

PairSPH_RHOSUM::PairSPH_RHOSUM(PDPS *ps) : Pair(ps)
{
	first = 1;
	newton_pair = 1;
	allocated = 0;
	cubic_flag = 0;
	quintic_flag = 0;
	h = 0.0;
	comm_forward = comm_reverse = 1;
	cut = NULL;
	cutsq = NULL;
	devSetflag = NULL;
	devCutsq = NULL;
}

/* ---------------------------------------------------------------------- */

PairSPH_RHOSUM::~PairSPH_RHOSUM()
{
	if (allocated) {
		memory->destroy(setflag);
		memory->destroy(cutsq);
		memory->destroy(rho0);
		memory->destroy(cut);

	}

}

/* ---------------------------------------------------------------------- */

void PairSPH_RHOSUM::allocate()
{
	allocated = 1;
	int n = particle->ntypes;

	memory->create(setflag, n + 1, n + 1, "pair:setflag");
	for (int i = 1; i <= n; i++)
	for (int j = i; j <= n; j++)
		setflag[i][j] = 0;
	memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
	memory->create(rho0, n + 1, "pair:rho0");
	memory->create(cut, n + 1, n + 1, "pair:cut");

	// pointer to transfer to GPU
	hostSetflag = (int *)malloc(TYPEMAX * TYPEMAX * sizeof(int));
	hostCutsq = (double *)malloc(TYPEMAX * TYPEMAX * sizeof(double));

}

/* ----------------------------------------------------------------------
Compute force for all paritcles
------------------------------------------------------------------------- */

void PairSPH_RHOSUM::compute(int eflag, int vflag)
{
	int i, j, ii, jj, inum, jnum, itype, jtype;
	double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

	int *ilist, *jlist, *numneigh, **firstneigh;
	double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, q;
	double rsq, rij_inv, tmp, wfd, delVdotDelR, mu, deltaE;
	MPI_Request request;
	MPI_Status status;
	if (eflag || vflag)
		ev_setup(eflag, vflag);
	//  else
	//    evflag = vflag_fdotr = 0;

	double **v = particle->vest;
	double **x = particle->x;
	double **f = particle->f;
	double *rho = particle->rho;
	double *mass = particle->mass;
	double *de = particle->de;
	double *drho = particle->drho;
	int *type = particle->type;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	double wf;
	//  int newton_pair = force->newton_pair;
	// check consistency of pair coefficients
	if (first) {
		for (i = 1; i <= particle->ntypes; i++) {
			for (j = 1; i <= particle->ntypes; i++) {
				if (cutsq[i][j] > 1.e-32) {
					if (!setflag[i][i] || !setflag[j][j]) {
						if (parallel->procid == 0) {
							printf(
								"SPH particle types %d and %d interact with cutoff=%g, but not all of their single particle properties are set.\n",
								i, j, sqrt(cutsq[i][j]));
						}
					}
				}
			}
		}
		first = 0;
	}

	inum = neighbor->neighlist->inum;
	ilist = neighbor->neighlist->ilist;
	numneigh = neighbor->neighlist->numneigh;
	firstneigh = neighbor->neighlist->firstneigh;

	// loop over neighbors of my particles
	if ((update->ntimestep % nstep) == 0) {

		//// initialize density with self-contribution,
		//for (i = 0; i < nlocal; i++) {
		//	itype = type[i];
		//	if (setflag[itype][itype]){
		//		imass = mass[itype];
		//		if (domain->dim == 3) {

		//			// Cubic spline kernel, 3d
		//			wf = a3D;
		//		}
		//		else {
		//			// Cubic spline kernel, 2d
		//			wf = a2D;
		//		}
		//		rho[i] = imass * wf;
		//	}
		//	
		//}
		////	set all ghost particle's rho zero to be computed
		//if (particle->nghost > 0){
		//	for (i = nlocal; i < nlocal + particle->nghost; i++){
		//		itype = type[i];
		//		if (setflag[itype][itype])
		//			rho[i] = 0.0;
		//	}
		//}
		//for (ii = 0; ii < inum; ii++) {
		//	i = ilist[ii];
		//	xtmp = x[i][0];
		//	ytmp = x[i][1];
		//	ztmp = x[i][2];
		//	itype = type[i];
		//	jlist = firstneigh[i];
		//	jnum = numneigh[i];

		//	for (jj = 0; jj < jnum; jj++) {
		//		j = jlist[jj];

		//		jtype = type[j];
		//		if (setflag[itype][itype] || setflag[jtype][jtype]){
		//			delx = xtmp - x[j][0];
		//			dely = ytmp - x[j][1];
		//			delz = ztmp - x[j][2];
		//			rsq = delx * delx + dely * dely + delz * delz;

		//			if (rsq < cutsq[itype][jtype]) {
		//				q = sqrt(rsq) / h;

		//				if (cubic_flag == 1){
		//					if (q < 1)
		//						wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
		//					else
		//						wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
		//				}
		//				else if (quintic_flag == 1)
		//					wf = (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (2 * q + 1);

		//				if (domain->dim == 3)
		//					wf = wf * a3D;
		//				else
		//					wf = wf * a2D;
		//				if (setflag[itype][itype]){
		//					rho[i] += mass[jtype] * wf;
		//				}
		//					
		//				if (setflag[jtype][jtype]){
		//					rho[j] += mass[itype] * wf;
		//				}


		//			}

		//		}

		//	}

		//}
		//// particles get its rho part in neighbor processors, so that they can get complete rho
		//parallel->reverse_comm_pair(this);
		////	after getting the correct rho, update its value to its neighbor processors
		//parallel->forward_comm_pair(this);

		//for (i = 0; i < 0; i++){
		//	itype = type[i];
		//	xtmp = x[i][0];
		//	ytmp = x[i][1];
		//	ztmp = x[i][2];	

		//	if (domain->dim == 3)
		//		drho[i] = a3D / rho[i];
		//	else
		//		drho[i] = a2D / rho[i];
		//	jlist = firstneigh[i];
		//	jnum = numneigh[i];
		//	for (jj = 0; jj < jnum; jj++) {
		//		j = jlist[jj];
		//		jtype = type[j];
		//		delx = xtmp - x[j][0];
		//		dely = ytmp - x[j][1];
		//		delz = ztmp - x[j][2];
		//		rsq = delx * delx + dely * dely + delz * delz;

		//		if (rsq < cutsq[itype][jtype]) {
		//			q = sqrt(rsq) / h;

		//			if (cubic_flag == 1){
		//				if (q < 1)
		//					wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
		//				else
		//					wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
		//			}
		//			else if (quintic_flag == 1)
		//				wf = (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (2 * q + 1);

		//			if (domain->dim == 3)
		//				wf = wf * a3D;
		//			else
		//				wf = wf * a2D;
		//		}

		//	}

		//}

		cudaError_t error_t;
		error_t = cudaMemcpy(particle->ptrHostRho, particle->devRho, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

		gpuComputerho << < int(nlocal + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >(particle->devCoordX, particle->devCoordY, particle->devCoordZ,
			neighbor->devPairtable, neighbor->devNumneigh, particle->devRho, particle->devMass,
			particle->devType, particle->devMask, h, nlocal, a3D, devSetflag, devCutsq);

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time, start, stop);
		time = time;
		//error_t = cudaMemcpy(hostSetflag, devSetflag, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
		//error_t = cudaMemcpy(hostCutsq, devCutsq, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
//		error_t = cudaMemcpy(particle->ptrHostRho, particle->devRho, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
		//error_t = cudaMemcpy(neighbor->hostNumneigh, neighbor->devNumneigh, particle->nlocal * sizeof(int), cudaMemcpyDeviceToHost);
		//error_t = cudaMemcpy(neighbor->hostPairtable, neighbor->devPairtable, particle->nlocal * NEIGHMAX * sizeof(int), cudaMemcpyDeviceToHost);
		//error_t = error_t;
	}

}



/* ----------------------------------------------------------------------
Setting for pair_style command
------------------------------------------------------------------------- */

void PairSPH_RHOSUM::set_style(int narg, char **arg)
{
	//	if (narg != 4)
	//		error->all(FLERR, "Illegal number of setting arguments for pair_style sph/idealgas");
	if (strcmp(arg[1], "Cubic") == 0)
		cubic_flag = 1;
	else if (strcmp(arg[1], "Quintic") == 0)
		quintic_flag = 1;
	else
		error->all(FLERR, "Wrong Kernel function");
	nstep = atoi(arg[2]);

}

/* ----------------------------------------------------------------------
Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairSPH_RHOSUM::set_coeff(int narg, char **arg)
{
	if (narg != 4)
		error->all(FLERR, "Incorrect args for pair_style sph/taitwater coefficients");
	if (!allocated)
		allocate();


	int ilo, ihi, jlo, jhi;
	force->bounds(arg[0], particle->ntypes, ilo, ihi);
	force->bounds(arg[1], particle->ntypes, jlo, jhi);

	double rho0_one = atof(arg[2]);
	double cut_one = atof(arg[3]);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		rho0[i] = rho0_one;
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			rho0[j] = rho0_one;
			cut[i][j] = 2 * cut_one;
			cutsq[i][j] = cut[i][j] * cut[i][j];
			setflag[i][j] = 1;
			count++;
			hostSetflag[i * TYPEMAX + j] = setflag[i][j];
			hostCutsq[i * TYPEMAX + j] = cutsq[i][j];
		}
	}

	h = cut_one;

	if (cubic_flag == 1){
		a2D = 10.0 / 7.0 / PI / h / h;
		a3D = 1.0 / PI / h / h / h;
	}
	else if (quintic_flag == 1){
		a2D = 7.0 / 4.0 / PI / h / h;
		a3D = 21.0 / 16.0 / PI / h / h / h;
	}

	if (count == 0)
		error->all(FLERR, "Incorrect args for pair coefficients");

	// setup for GPU parameters
	cudaError_t cudaStatus;
	cudaMalloc(&devSetflag, TYPEMAX * TYPEMAX * sizeof(int));
	cudaMalloc(&devCutsq, TYPEMAX * TYPEMAX * sizeof(double));
	cudaMemcpy(devSetflag, hostSetflag, TYPEMAX * TYPEMAX * sizeof(int), cudaMemcpyHostToDevice);
	cudaStatus = cudaMemcpy(devCutsq, hostCutsq, TYPEMAX * TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
}

/* ----------------------------------------------------------------------
init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairSPH_RHOSUM::init_one(int i, int j) {

	if (setflag[i][j] == 0) {
		error->all(FLERR, "Not all pair sph/taitwater coeffs are set");
	}
	force->type2pair[i][j] = pair_id;

	cut[j][i] = cut[i][j];
	cutsq[j][i] = cutsq[i][j];
	setflag[j][i] = setflag[i][j];
	force->type2pair[j][i] = force->type2pair[i][j];


}

/* ---------------------------------------------------------------------- */

double PairSPH_RHOSUM::single(int i, int j, int itype, int jtype,
	double rsq, double factor_coul, double factor_lj, double &fforce) {
	fforce = 0.0;
	return 0.0;
}

/* ----------------------------------------------------------------------
	pack particle's rho to neighbor processors during reverse communication
------------------------------------------------------------------------- */

int PairSPH_RHOSUM::pack_reverse_comm(int n, int first, double *buf) {
	int i, m, last;
	double *rho = particle->rho;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		buf[m++] = rho[i];
	}
	return m;
}

/* ----------------------------------------------------------------------
	unpack particle's rho from neighbor processors during reverse communication
------------------------------------------------------------------------- */

void PairSPH_RHOSUM::unpack_reverse_comm(int n, int *list, double *buf) {
	int i, m, j;
	double *rho = particle->rho;
	int *type = particle->type;
	int jtype;
	m = 0;
	
	for (i = 0; i < n; i++) {
		j = list[i];
		jtype = type[j];
		if (setflag[jtype][jtype])
			rho[j] += buf[m++];
		else
			m++;
	}
}

/* ----------------------------------------------------------------------
	pack particle's rho to neighbor processors during forward communication
------------------------------------------------------------------------- */

int PairSPH_RHOSUM::pack_forward_comm(int n, int *list,  double *buf) {
	int i, m, j; 
	double *rho = particle->rho;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		buf[m++] = rho[j];
	}
	return m;
}

/* ----------------------------------------------------------------------
	unpack particle's rho from neighbor processors during forward communication
------------------------------------------------------------------------- */

void PairSPH_RHOSUM::unpack_forward_comm(int n, int first, double *buf) {
	int i, m, last;
	double *rho = particle->rho;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		rho[i] = buf[m++];
	}

}

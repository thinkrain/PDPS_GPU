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
#include "pair_sph_taitwater.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"
#include "group.h"
#include "mpi.h"
#include "timer.h"

#include "pdps_cuda.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
using namespace PDPS_NS;


#define DELTA 1
#define EPSILON 1.0e-10
#define PI 3.1416

enum{ ARTIFICIAL, LAMINAR, SPS };

__global__ void gpuComputesphforce(double *devCoordX, double *devCoordY, double *devCoordZ, int *devPairtable, int *devNumneigh,
	double *devRho,  double *devMass, int *devType, int *devMask, const double h, const int nlocal, const double a3D,  int *devSetflag, 
	double *devCutsq, double *devRho0, double *devB, const int visc_flag, double *devVeloX, double *devVeloY, double *devVeloZ, 
	double *devVisc, double *devSoundspeed,double *devForceX, double *devForceY, double *devForceZ){
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	float wf, xtemp, ytemp, ztemp, rsq, delx, dely, delz, q, tmp, fi, fj, irho, jrho, rij_inv, delVdotDelR;
	float vxtmp, vytmp, vztmp, mu, fviscx, fviscy, fviscz, wfd, fpair, imass, jmass;
	unsigned int j, jj, itype, jtype, jnum;
	__shared__ float mass[TYPEMAX];
	__shared__ float rho0[TYPEMAX];
	__shared__ float B[TYPEMAX];
	__shared__ int setflag[TYPEMAX * TYPEMAX];
	__shared__ float cutsq[TYPEMAX * TYPEMAX];
	__shared__ float viscosity[TYPEMAX * TYPEMAX];
	__shared__ float soundspeed[TYPEMAX];
	for (int tid = 0; tid < TYPEMAX; tid++){
		mass[tid] = devMass[tid];
		rho0[tid] = devRho0[tid];
		B[tid] = devB[tid];
		soundspeed[tid] = devSoundspeed[tid];
		for (j = 0; j < TYPEMAX; j++){
			setflag[tid * TYPEMAX + j] = devSetflag[tid * TYPEMAX + j];
			cutsq[tid * TYPEMAX + j] = devCutsq[tid * TYPEMAX + j];
			viscosity[tid * TYPEMAX + j] = devVisc[tid * TYPEMAX + j];
		}

	}

	
	for (i = i; i < nlocal; i += blockDim.x * gridDim.x){

	
			itype = devType[i];
			if (setflag[itype * TYPEMAX + itype] == 1){
				irho = devRho[i];
				xtemp = devCoordX[i];
				ytemp = devCoordY[i];
				ztemp = devCoordZ[i];
				vxtmp = devVeloX[i];
				vytmp = devVeloY[i];
				vztmp = devVeloZ[i];
				imass = mass[itype];
				jnum = devNumneigh[i];
				tmp = irho / rho0[itype];
				fi = tmp * tmp * tmp;
				fi = B[itype] * (fi * fi *tmp - 1.0) / (irho * irho);
				for (jj = 0; jj < jnum; jj++){
					j = devPairtable[i * NEIGHMAX + jj];
					jtype = devType[j];
					if (setflag[itype * TYPEMAX + jtype]){
						delx = xtemp - devCoordX[j];
						dely = ytemp - devCoordY[j];
						delz = ztemp - devCoordZ[j];
						rsq = delx * delx + dely * dely + delz * delz;
						jmass = mass[jtype];
						jrho = devRho[j];
						if (rsq < cutsq[itype * TYPEMAX + jtype]){
							q = sqrt(rsq) / h;
							rij_inv = 1.0 / sqrt(rsq);
							//	Cubic Spline
							if (q < 1)
								wfd = -3 * q + 2.25 * q * q;
							else
								wfd = -0.75 * (2 - q) * (2 - q);
							wfd = wfd * a3D / h;

							tmp = jrho / rho0[jtype];
							fj = tmp * tmp * tmp;
							fj = B[jtype] * (fj * fj * tmp - 1.0) / (jrho * jrho);
							delVdotDelR = delx * (vxtmp - devVeloX[j]) + dely * (vytmp - devVeloY[j])
								+ delz * (vztmp - devVeloZ[j]);

							if (visc_flag == 0){
								// artificial viscosity (Monaghan 1992)
								if (delVdotDelR < 0) {
									mu = h * delVdotDelR / (rsq + 0.01 * h * h);
									fviscx = fviscy = fviscz = -viscosity[itype * TYPEMAX + jtype] * (soundspeed[itype]
										+ soundspeed[jtype]) * mu / (irho + jrho);
									fviscx *= -imass * jmass * wfd * delx * rij_inv;
									fviscy *= -imass * jmass * wfd * dely * rij_inv;
									fviscz *= -imass * jmass * wfd * delz * rij_inv;

								}
								else {
									fviscx = fviscy = fviscz = 0.0;
								}
							}
							else if (visc_flag == 1){
								// Laminar viscosity (Lo and Shao 2002)
								mu = 2 * viscosity[itype * TYPEMAX + jtype] / (irho + jrho) / (rsq + 0.01 * h * h) * wfd / rij_inv;
								fviscx = mu * imass * jmass * (vxtmp - devVeloX[j]);
								fviscy = mu * imass * jmass * (vytmp - devVeloY[j]);
								fviscz = mu * imass * jmass * (vztmp - devVeloZ[j]);

							}
					
							// total pair force & thermal energy increment
							fpair = -imass * jmass * (fi + fj) * wfd;

							devForceX[i] += delx * fpair * rij_inv + fviscx;
							devForceY[i] += dely * fpair * rij_inv + fviscy;
							devForceZ[i] += delz * fpair * rij_inv + fviscz;



						}		//  rsq < cutsq[itype][jtype]

					}	// setflag[itype * 10 + jtype]
				}	// j < jnum
			}	// setflag[itype * TYPEMAX + itype]

	}


}

/* ---------------------------------------------------------------------- */

PairSPH_TAITWATER::PairSPH_TAITWATER(PDPS *ps) : Pair(ps)
{

	 first = 1;
	 newton_pair = 1;
	 allocated = 0;
	 cubic_flag = 0;
	 quintic_flag = 0;
	 h = 0.0;
	 cut = NULL;
	 cutsq = NULL;
}

/* ---------------------------------------------------------------------- */

PairSPH_TAITWATER::~PairSPH_TAITWATER()
{
	if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(rho0);
    memory->destroy(soundspeed);
    memory->destroy(B);
    memory->destroy(viscosity);
  }

}

/* ---------------------------------------------------------------------- */

void PairSPH_TAITWATER::allocate()
{
	 allocated = 1;
  int n = particle->ntypes;

  // pointer to transfer to GPU
  hostSetflag = (int *)malloc(TYPEMAX * TYPEMAX * sizeof(int));
  hostCutsq = (double *)malloc(TYPEMAX * TYPEMAX * sizeof(double));
  hostVisc = (double *)malloc(TYPEMAX * TYPEMAX * sizeof(double));

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
  for (int j = i; j <= n; j++){
	  setflag[i][j] = 0;
	  hostSetflag[i * TYPEMAX + j] = 0;
	  hostCutsq[i * TYPEMAX + j] = 0;
	  hostVisc[i * TYPEMAX + j] = 0;
  }


  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(rho0, n + 1, "pair:rho0");
  memory->create(soundspeed, n + 1, "pair:soundspeed");
  memory->create(B, n + 1, "pair:B");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");


 

}

/* ----------------------------------------------------------------------
   Compute force for all paritcles
------------------------------------------------------------------------- */

void PairSPH_TAITWATER::compute(int eflag, int vflag)
{


	int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fviscx, fviscy, fviscz, q;
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
 // for (ii = 0; ii < inum; ii++) {
 //   i = ilist[ii];
 //   xtmp = x[i][0];
 //   ytmp = x[i][1];
 //   ztmp = x[i][2];
 //   vxtmp = v[i][0];
 //   vytmp = v[i][1];
 //   vztmp = v[i][2];
 //   itype = type[i];
 //   jlist = firstneigh[i];
 //   jnum = numneigh[i];

 //   imass = mass[itype];
 //   // compute pressure of particle i with Tait EOS
 //   tmp = rho[i] / rho0[itype];
 //   fi = tmp * tmp * tmp;
	//fi = B[itype] * (fi * fi * tmp - 1.0) / (rho[i] * rho[i]);

 //   for (jj = 0; jj < jnum; jj++) {
 //     j = jlist[jj];
	//  jtype = type[j];

	//  if (setflag[itype][jtype]){
	//	  delx = xtmp - x[j][0];
	//	  dely = ytmp - x[j][1];
	//	  delz = ztmp - x[j][2];
	//	  rsq = delx * delx + dely * dely + delz * delz;
	//	  jmass = mass[jtype];

	//	  if (rsq < cutsq[itype][jtype]){
	//		  q = sqrt(rsq) / h;
	//		  rij_inv = 1.0 / sqrt(rsq);
	//		  if (cubic_flag == 1){
	//			  if (q < 1)
	//				  wfd = -3 * q + 2.25 * q * q;
	//			  else
	//				  wfd = -0.75 * (2 - q) * (2 - q);
	//		  }
	//		  else if (quintic_flag == 1)
	//			  wfd = -2 * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (2 * q + 1) + 2 * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0);

	//		  if (domain->dim == 3)
	//			  wfd = wfd * a3D / h;
	//		  else
	//			  wfd = wfd * a2D / h;

	//		  // compute pressure  of particle j with Tait EOS
	//		  tmp = rho[j] / rho0[jtype];
	//		  fj = tmp * tmp * tmp;
	//		  fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[j] * rho[j]);
	//		  // dot product of velocity delta and distance vector
	//		  delVdotDelR = delx * (vxtmp - v[j][0]) + dely * (vytmp - v[j][1])
	//			  + delz * (vztmp - v[j][2]);

	//		  if (visc_flag == ARTIFICIAL){
	//			  // artificial viscosity (Monaghan 1992)
	//			  if (delVdotDelR < 0.) {
	//				  mu = h * delVdotDelR / (rsq + 0.01 * h * h);
	//				  fviscx = fviscy = fviscz = -viscosity[itype][jtype] * (soundspeed[itype]
	//					  + soundspeed[jtype]) * mu / (rho[i] + rho[j]);
	//				  fviscx *= -imass * jmass * wfd * delx / rij_inv;
	//				  fviscy *= -imass * jmass * wfd * dely / rij_inv;
	//				  fviscz *= -imass * jmass * wfd * delz / rij_inv;

	//			  }
	//			  else {
	//				  fviscx = fviscy = fviscz = 0.;
	//			  }
	//		  }
	//		  else if (visc_flag == LAMINAR){
	//			  // Laminar viscosity (Lo and Shao 2002)
	//			  mu = 2 * viscosity[itype][jtype] / (rho[i] + rho[j]) / (rsq + 0.01 * h * h) * wfd / rij_inv;
	//			  fviscx = mu * imass * jmass * (vxtmp - v[j][0]);
	//			  fviscy = mu * imass * jmass * (vytmp - v[j][1]);
	//			  fviscz = mu * imass * jmass * (vztmp - v[j][2]);

	//		  }


	//		  // total pair force & thermal energy increment
	//		  fpair = -imass * jmass * (fi + fj) * wfd;
	//		  deltaE = -0.5 * fpair * delVdotDelR;
	//		  f[i][0] += delx * fpair * rij_inv + fviscx;
	//		  f[i][1] += dely * fpair * rij_inv + fviscy;
	//		  f[i][2] += delz * fpair * rij_inv + fviscz;
	//		  // and change in density
	//		  drho[i] += jmass * delVdotDelR * wfd;
	//		  // change in thermal energy
	//		  de[i] += deltaE;

	//		  if (newton_pair || j < nlocal) {
	//			  f[j][0] -= delx * fpair * rij_inv + fviscx;
	//			  f[j][1] -= dely * fpair * rij_inv + fviscy;
	//			  f[j][2] -= delz * fpair * rij_inv + fviscz;
	//			  de[j] += deltaE;
	//			  drho[j] += imass * delVdotDelR * wfd;
	//		  }

	//		  if (evflag)
	//			  ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
	//	  }

	//  } // setflag[itype][jtype]

 //    
 //   }

 // }
 cudaError_t error_t;
  //force->clear();
  /*error_t = cudaMemcpy(neighbor->hostForceX, particle->devForceX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
  error_t = cudaMemcpy(neighbor->hostForceZ, particle->devForceZ, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
  error_t = cudaMemcpy(particle->ptrHostRho, particle->devRho, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);*/
  //error_t = cudaMemcpy(neighbor->hostForceX, particle->devForceX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
  //error_t = cudaMemcpy(neighbor->hostForceY, particle->devForceY, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
  //error_t = cudaMemcpy(neighbor->hostForceZ, particle->devForceZ, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	 cudaEvent_t start, stop;
	 float time;
	 cudaEventCreate(&start);
	 cudaEventCreate(&stop);
	 cudaEventRecord(start, 0);

	 gpuComputesphforce << < int(nlocal + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >(particle->devCoordX, particle->devCoordY, particle->devCoordZ,
	  neighbor->devPairtable, neighbor->devNumneigh, particle->devRho, particle->devMass, particle->devType,
	  particle->devMask, h, nlocal, a3D, devSetflag, devCutsq, devRho0, devB, visc_flag,
	  particle->devVestX, particle->devVestY, particle->devVestZ, devVisc, devSoundspeed,
	  particle->devForceX, particle->devForceY, particle->devForceZ);

	  cudaEventRecord(stop, 0);
	  cudaEventSynchronize(stop);
	  cudaEventElapsedTime(&time, start, stop);
 /* error_t = cudaMemcpy(neighbor->hostForceX, particle->devForceX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
  error_t = cudaMemcpy(neighbor->hostForceY, particle->devForceY, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
  error_t = cudaMemcpy(neighbor->hostForceZ, particle->devForceZ, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);*/
  //error_t = cudaMemcpy(hostCutsq, devCutsq, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
 /* error_t = cudaMemcpy(neighbor->hostNumneigh, neighbor->devNumneigh, particle->nlocal * sizeof(int), cudaMemcpyDeviceToHost);
  error_t = cudaMemcpy(neighbor->hostPairtable, neighbor->devPairtable, particle->nlocal * NEIGHMAX * sizeof(int), cudaMemcpyDeviceToHost);
  error_t = cudaMemcpy(hostSetflag, devSetflag, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
  error_t = cudaMemcpy(hostCutsq, devSetflag, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
  error_t = cudaMemcpy(particle->ptrHostRho, particle->devRho, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
  error_t = cudaMemcpy(neighbor->hostForceX, particle->devForceX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
  error_t = cudaMemcpy(neighbor->hostForceZ, particle->devForceZ, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);*/

 // //particle->TransferG2C();
 
}



/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairSPH_TAITWATER::set_style(int narg, char **arg)
{
//	if (narg != 4)
//		error->all(FLERR, "Illegal number of setting arguments for pair_style sph/idealgas");
	if (strcmp(arg[1], "Cubic") == 0)
		cubic_flag = 1;  
	else if (strcmp(arg[1], "Quintic") == 0)
		quintic_flag = 1;
	else
		error->all(FLERR, "Wrong Kernel function");
	//	viscosity type
	if (strcmp(arg[2], "Artificial") == 0)
		visc_flag = ARTIFICIAL;
	else if (strcmp(arg[2], "Laminar") == 0)
		visc_flag = LAMINAR;
	else if (strcmp(arg[2], "SPS") == 0)
		visc_flag = SPS;
	else
		error->all(FLERR, "Wrong Viscosity function");

}

/* ----------------------------------------------------------------------
                              Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairSPH_TAITWATER::set_coeff(int narg, char **arg)
{
   if (narg != 6)
    error->all(FLERR, "Incorrect args for pair_style sph/taitwater coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;

  force->bounds(arg[0], particle->ntypes, ilo, ihi);
  force->bounds(arg[1], particle->ntypes, jlo, jhi);

  double rho0_one = atof(arg[2]);
  double soundspeed_one = atof(arg[3]);
  double viscosity_one = atof(arg[4]);
  double cut_one = atof(arg[5]);
  double B_one = soundspeed_one * soundspeed_one * rho0_one / 7.0;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rho0[i] = rho0_one;
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    for (int j = MAX(jlo,i); j <= jhi; j++) {
		rho0[j] = rho0_one;
		soundspeed[j] = soundspeed_one;
		B[j] = B_one;
		viscosity[i][j] = viscosity_one;
		cut[i][j] = 2 * cut_one;
		cutsq[i][j] = cut[i][j] * cut[i][j];
		setflag[i][j] = 1;

		hostSetflag[i * TYPEMAX + j] = setflag[i][j];
		hostCutsq[i * TYPEMAX + j] = cutsq[i][j];
		hostVisc[i * TYPEMAX + j] = viscosity[i][j];
      count++;
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
    error->all(FLERR,"Incorrect args for pair coefficients");

  // setup for GPU parameters
  cudaError_t cudaStatus;

  cudaMalloc(&devVisc, TYPEMAX * TYPEMAX * sizeof(double));
  cudaMalloc(&devB, TYPEMAX * sizeof(double));
  cudaMalloc(&devSoundspeed, TYPEMAX * sizeof(double));
  cudaMalloc(&devRho0, TYPEMAX * sizeof(double));
  cudaMalloc(&devSetflag, TYPEMAX * TYPEMAX * sizeof(int));
  cudaMalloc(&devCutsq, TYPEMAX * TYPEMAX * sizeof(double));

  cudaMemcpy(devVisc, hostVisc, TYPEMAX * TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(devB, B, TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(devSoundspeed, soundspeed, TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(devRho0, rho0, TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(devSetflag, hostSetflag, TYPEMAX * TYPEMAX * sizeof(int), cudaMemcpyHostToDevice);
  cudaStatus = cudaMemcpy(devCutsq, hostCutsq, TYPEMAX * TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairSPH_TAITWATER::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Not all pair sph/taitwater coeffs are set");
  }
  force->type2pair[i][j] = pair_id;

  cut[j][i] = cut[i][j];
  cutsq[j][i] = cutsq[i][j];
  setflag[j][i] = setflag[i][j];
  viscosity[j][i] = viscosity[i][j];
  force->type2pair[j][i] = force->type2pair[i][j];

}

/* ---------------------------------------------------------------------- */

double PairSPH_TAITWATER::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;
  return 0.0;
}

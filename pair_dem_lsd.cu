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
#include <iostream>
#include <fstream>
#include <iomanip>

#include "error.h"
#include "domain.h"
#include "force.h"
#include "pair_list.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair_dem_lsd.h"
#include "parallel.h"
#include "particle.h"
#include "psmath.h"
#include "update.h"

#include "output.h"

using namespace PDPS_NS;
using namespace PsMath_NS;

#define DELTA 1
#define EPSILON 1.0e-10
#define float_double double 

__global__ void gpuComputedemforce(
    double *devCoordX, double *devCoordY, double *devCoordZ,
    double *devVeloX, double *devVeloY, double *devVeloZ,  
	int *devPairtable, int *devNumneigh,
//	int *devContactFlag, double *devDrijtx, double *devDrijty, double *devDrijtz,
	double *devMass, double *devRadius, int *devType, int *devMask, const int nlocal, const double dt, 
	const int radius_flag, const int torque_flag, int *devSetflag, 
	double *devCutsq, double *devKn, double *devKt, double *devCn, double *devCt, double *devMu,
//	double *devOmegaX, double *devOmegaY, double *devOmegaZ,
//	double *devTorqueX, double *devTorqueY, double *devTorqueZ,
	double *devForceX, double *devForceY, double *devForceZ){
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	float wf, ix, iy, iz, rsq, rij, radius_cut, rijx, rijy, rijz, q, tmp, fi, fj, irho, jrho, rij_inv, delVdotDelR;
	float ivx, ivy, ivz, fviscx, fviscy, fviscz, wfd, fpair, imass, jmass, nx, ny, nz;
	float drijn, drijnx, drijny, drijnz, vijx, vijy, vijz, Li, Lj;
	float omegainijx, omegainijy, omegainijz, omegajnijx, omegajnijy, omegajnijz;
	float vijn, vijnx, vijny, vijnz, vijt, vijtx, vijty, vijtz, vijt_inv, tx, ty, tz;
	float drijtx, drijty, drijtz, temp, torqueijx, torqueijy, torqueijz;
	float fijn, fijnx, fijny, fijnz, fijt, fijtx, fijty, fijtz;
	unsigned int j, jj, itype, jtype, jnum;
	int contact_flag;
	
	__shared__ float mass[TYPEMAX];
	__shared__ int setflag[TYPEMAX * TYPEMAX];
	__shared__ float cutsq[TYPEMAX * TYPEMAX];
	__shared__ float kn[TYPEMAX * TYPEMAX];
	__shared__ float kt[TYPEMAX * TYPEMAX];
	__shared__ float Cn[TYPEMAX * TYPEMAX];
	__shared__ float Ct[TYPEMAX * TYPEMAX];
	__shared__ float mu[TYPEMAX * TYPEMAX];
	
	for (int tid = 0; tid < TYPEMAX; tid++){
		mass[tid] = devMass[tid];
		for (j = 0; j < TYPEMAX; j++){
			setflag[tid * TYPEMAX + j] = devSetflag[tid * TYPEMAX + j];
			cutsq[tid * TYPEMAX + j] = devCutsq[tid * TYPEMAX + j];
			kn[tid * TYPEMAX + j] = devKn[tid * TYPEMAX + j];
			kt[tid * TYPEMAX + j] = devKt[tid * TYPEMAX + j];
			Cn[tid * TYPEMAX + j] = devCn[tid * TYPEMAX + j];
			Ct[tid * TYPEMAX + j] = devCt[tid * TYPEMAX + j];
			mu[tid * TYPEMAX + j] = devMu[tid * TYPEMAX + j];
		}

	}
	
	for (i = i; i < nlocal; i += blockDim.x * gridDim.x){
	
		itype = devType[i];
//		if (setflag[itype * TYPEMAX + itype] == 1){
			ix = devCoordX[i];
			iy = devCoordY[i];
			iz = devCoordZ[i];
			ivx = devVeloX[i];
			ivy = devVeloY[i];
			ivz = devVeloZ[i];
			//imass = mass[itype];
			jnum = devNumneigh[i];
			for (jj = 0; jj < jnum; jj++){
				j = devPairtable[i * NEIGHMAX + jj];
				jtype = devType[j];
				if (setflag[itype * TYPEMAX + jtype]){
					rijx = ix - devCoordX[j];
					rijy = iy - devCoordY[j];
					rijz = iz - devCoordZ[j];
					rsq = rijx * rijx + rijy * rijy + rijz * rijz;
					//jmass = mass[jtype];
					if (rsq < cutsq[itype * TYPEMAX + jtype]){
						rij = sqrt(rsq);
						if (rij < EPSILON) continue;    // rij cannot be 0.0 in DEM systems
						if (radius_flag == 1) {
							radius_cut = devRadius[i] + devRadius[j];
							//if (rij >= radius_cut) {
							//	devContactFlag[i * NEIGHMAX + jj] = 0;
							//	devDrijtx[i * NEIGHMAX + jj] = 0.0;
							//	devDrijty[i * NEIGHMAX + jj] = 0.0;
							//	devDrijtz[i * NEIGHMAX + jj] = 0.0;
							//	continue;
							//}
						}
						
						rij_inv = 1.0/rij;
						// nij = rij/|rij|
						nx = rijx * rij_inv;              
						ny = rijy * rij_inv;
						nz = rijz * rij_inv;

						// overlap between two particles = rcut - (rij . nij) = rcut - |rij|
						// may need to be changed in the future by replaicng with particles's radius
						if (radius_flag) drijn = radius_cut - rij;
						else drijn = sqrt(cutsq[itype * TYPEMAX + jtype]) - rij;
						drijnx = drijn * nx;
						drijny = drijn * ny;
						drijnz = drijn * nz;
						
						vijx = ivx - devVeloX[j];
						vijy = ivy - devVeloY[j];
						vijz = ivz - devVeloZ[j];
						
//						if (radius_flag) {
//							// distance from the center of particle to the contact point
//							Li = (rij*rij + devRadius[i] * devRadius[i] - devRadius[j] * devRadius[j]) / (2 * rij);
//							Lj = rij - Li;
//							// vij = vi - vj + (wi cross Ri) - (wj cross Rj) (relative velocity)  
//							omegainijx = devOmegaY[i]*nz - devOmegaZ[i]*ny;
//							omegainijy = devOmegaZ[i]*nx - devOmegaX[i]*nz;
//							omegainijz = devOmegaX[i]*ny - devOmegaY[i]*nx;
//							omegajnijx = devOmegaY[j]*nz - devOmegaZ[j]*ny;
//							omegajnijy = devOmegaZ[j]*nx - devOmegaX[j]*nz;
//							omegajnijz = devOmegaX[j]*ny - devOmegaY[j]*nx;
//	
//							vijx += -Li*omegainijx - Lj*omegajnijx;
//							vijy += -Li*omegainijy - Lj*omegajnijy;
//							vijz += -Li*omegainijz - Lj*omegajnijz;
//						}

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
						/*contact_flag = devContactFlag[i * NEIGHMAX + jj];
						if (contact_flag == 0) {*/
							if (fabs(vijn) < EPSILON) temp = dt;
							else temp = MIN(fabs(drijn/vijn), dt);
							drijtx = vijtx * temp;
							drijty = vijty * temp;
							drijtz = vijtz * temp;
						//	devContactFlag[i * NEIGHMAX + jj] = 1;
						//}
						// during the same contact
						//else {
						//	drijtx = devDrijtx[i * NEIGHMAX + jj];
						//	drijty = devDrijty[i * NEIGHMAX + jj];
						//	drijtz = devDrijtz[i * NEIGHMAX + jj];
						//	
						//	// update the tang. disp. for the next time step
						//	drijtx = drijtx + vijtx*dt;
						//	drijty = drijty + vijty*dt;
						//	drijtz = drijtz + vijtz*dt;
	
						//	temp = drijtx*nx + drijty*ny + drijtz*nz;
	
						//	drijtx = drijtx - temp*nx;
						//	drijty = drijty - temp*ny;
						//	drijtz = drijtz - temp*nz;
						//}
	
						//devDrijtx[i * NEIGHMAX + jj] = drijtx;
						//devDrijty[i * NEIGHMAX + jj] = drijty;
						//devDrijtz[i * NEIGHMAX + jj] = drijtz;

						// fijn = -kn*delta_rijn - Cn*vijn
						int itype_jtype = itype * TYPEMAX + jtype;
						fijnx = kn[itype_jtype]*drijnx - Cn[itype_jtype]*vijnx;
						fijny = kn[itype_jtype]*drijny - Cn[itype_jtype]*vijny;
						fijnz = kn[itype_jtype]*drijnz - Cn[itype_jtype]*vijnz;
	
						fijn = sqrt(fijnx*fijnx + fijny*fijny + fijnz*fijnz);
						
						// tangential force
						// fijt = -kt*delta_rijt - Ct*vijt
						fijtx = -kt[itype_jtype] * drijtx - Ct[itype_jtype] * vijtx;
						fijty = -kt[itype_jtype] * drijty - Ct[itype_jtype] * vijty;
						fijtz = -kt[itype_jtype] * drijtz - Ct[itype_jtype] * vijtz;
	
						fijt = sqrt(fijtx*fijtx + fijty*fijty + fijtz*fijtz); 
						temp = fabs(mu[itype_jtype]*fijn);
						// |fijt| > mu*|fnij|
						if (fijt > temp) {
							fijtx = -temp * tx;
							fijty = -temp * ty;
							fijtz = -temp * tz;
						}

						devForceX[i] += (fijnx + fijtx);
						devForceY[i] += (fijny + fijty);
						devForceZ[i] += (fijnz + fijtz);

//						if (torque_flag) {
//							torqueijx = ny*fijtz - nz*fijty;
//							torqueijy = nz*fijtx - nx*fijtz;
//							torqueijz = nx*fijty - ny*fijtx;
//	
//							devTorqueX[i] += -Li * torqueijx;
//							devTorqueY[i] += -Li * torqueijy;
//							devTorqueZ[i] += -Li * torqueijz;
//						}

					} // if (rij < cut[itype][jtype])
					//else {
					//	devContactFlag[i * NEIGHMAX + jj] = 0;
					//	devDrijtx[i * NEIGHMAX + jj] = 0.0;
					//	devDrijty[i * NEIGHMAX + jj] = 0.0;
					//	devDrijtz[i * NEIGHMAX + jj] = 0.0;
					//} // if (rij >= cut[itype][jtype])
						
				}	// setflag[itype * 10 + jtype]
			}	// j < jnum
		//}	// setflag[itype * TYPEMAX + itype]
		
	}

}

/* ---------------------------------------------------------------------- */

PairDEMLsd::PairDEMLsd(PDPS *ps) : Pair(ps)
{
	kn = Cn = kt = Ct = mu = NULL;
	pair_list = NULL;
	
	allocated = 0;

	tbsize = neighbor->pgsize;

	newton_pair = 1;

	pair_list = new PairList(ps,tbsize);        // initialize hashtable
	pair_list->ifields = 1;
	pair_list->dfields = 3;
	pair_list->pair_id = pair_id;

	rot_flag = 1;                               // By default, rotation is considered
//	if (particle->atomic_flag == 1) error->all(FLERR, "Illegal particle type");
}

/* ---------------------------------------------------------------------- */

PairDEMLsd::~PairDEMLsd()
{
	if (allocated) {
		memory->destroy(cut);
		memory->destroy(cutsq);
		memory->destroy(kn);
		memory->destroy(Cn);
		memory->destroy(kt);
		memory->destroy(Ct);
		memory->destroy(mu);
		memory->destroy(setflag);
	}

	delete pair_list;
	pair_list = NULL;
}

/* ----------------------------------------------------------------------
                              Setup
------------------------------------------------------------------------- */

void PairDEMLsd::setup()
{
	//int size = neighbor->neighlist->maxpage * neighbor->neighlist->pgsize;
	int size = neighbor->neighlist->pgsize;
	pair_list->init_hash(size);
}

/* ----------------------------------------------------------------------
                              Allocate
------------------------------------------------------------------------- */

void PairDEMLsd::allocate()
{
	allocated = 1;

	int n = particle->ntypes;

	// pointer to transfer to GPU
    hostSetflag = (int *)malloc(TYPEMAX * TYPEMAX * sizeof(int));
    hostCutsq = (double *)malloc(TYPEMAX * TYPEMAX * sizeof(double));
    hostKn = (double *)malloc(TYPEMAX * TYPEMAX * sizeof(double));
	hostKt = (double *)malloc(TYPEMAX * TYPEMAX * sizeof(double));
	hostCn = (double *)malloc(TYPEMAX * TYPEMAX * sizeof(double));
	hostCt = (double *)malloc(TYPEMAX * TYPEMAX * sizeof(double));
	hostMu = (double *)malloc(TYPEMAX * TYPEMAX * sizeof(double));
	
	setflag = memory->create(setflag, n+1, n+1, "PairDEM: setflag");
    for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++){
  	    setflag[i][j] = 0;
  	    hostSetflag[i * TYPEMAX + j] = 0;
  	    hostCutsq[i * TYPEMAX + j] = 0;
  	    hostKn[i * TYPEMAX + j] = 0;
		hostKt[i * TYPEMAX + j] = 0;
		hostCn[i * TYPEMAX + j] = 0;
		hostCt[i * TYPEMAX + j] = 0;
		hostMu[i * TYPEMAX + j] = 0;
    }
  
	cut = memory->create(cut, n+1, n+1, "PairDEM: cut");
	cutsq = memory->create(cutsq, n+1, n+1, "PairDEM: cutsq");
	kn = memory->create(kn, n+1, n+1, "PairDEM: kn");
	Cn = memory->create(Cn, n+1, n+1, "PairDEM: Cn");
	kt = memory->create(kt, n+1, n+1, "PairDEM: kt");
	Ct = memory->create(Ct, n+1, n+1, "PairDEM: Ct");
	mu = memory->create(mu, n+1, n+1, "PairDEM: mu");
	
	 // initialize setflag for pair i and j
	for (int i = 0; i <= n; i++)
	for (int j = 0; j <= n; j++) {
		setflag[i][j] = 0;
	}
}

/* ----------------------------------------------------------------------
   Compute force for all paritcles
------------------------------------------------------------------------- */

void PairDEMLsd::compute(int eflag, int vflag)
{
	int i, j, ii, jj, inum, jnum;        
	int itag, jtag, itype, jtype;
	double ix, iy, iz, ivx, ivy, ivz;                    // pos and vel for particle i
	double rij, rijsq, rij_inv, rijx, rijy, rijz;        // relative position
	double drijn, drijnx, drijny, drijnz;                // normal displacement
	double drijtx, drijty, drijtz;                       // tangential displacement
	double vijx, vijy, vijz;                             // relative velocity: vij = vj - vi
	double vijn, vijnx, vijny, vijnz;                    // relative velocity along normal direction
	double vijt, vijt_inv, vijtx, vijty, vijtz;          // relative velocity along tangential direction
	double fijn, fijnx, fijny, fijnz;                    // normal force
	double fijt_a[3], fijt, fijtx, fijty, fijtz;         // tangential force	
	double n[3], nx, ny, nz;                             // unit vector along normal direction
	double tx, ty, tz;                                   // unit vector along tangential direction
	double evdwl;
	double omegainij[3], omegajnij[3];                   // omega_i cross nij
	double torqueij[3];                                  // torque
	double Li, Lj;                                       // distance from the center of particle to the contact point             

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double *radius = particle->radius;
	double **omega = particle->omega;
	double **torque = particle->torque;
	int *type = particle->type;
	int *tag = particle->tag;
	int radius_flag = particle->radius_flag;
	int nlocal = particle->nlocal;
	int nall = nlocal + particle->nghost;               // number of particles stored in this processor 
	                                                    // including ghost particles
	double dt = update->dt;

	// neighbor class related
	int *ilist;                                    
	int *jlist;
	int **list;
	int **firstneigh;
	int *numneigh;

	ilist = neighbor->neighlist->ilist;
	inum = neighbor->neighlist->inum;
	firstneigh = neighbor->neighlist->firstneigh;
	numneigh = neighbor->neighlist->numneigh;

	evdwl = 0.0;
	if (eflag || vflag) ev_setup(eflag,vflag);

	// if the neighbor has been rebuilt, the pair_list needs to be rebuilt
//	int nflag = neighbor->nflag;
//	if (nflag == 1) {
//		set_pairlist();
//	}

//	PairList::HashPair **hash_pair = pair_list->hash_pair;
//	int contact_flag;
//	
//	char str[BYTES];
//	// loop for all local owned particles 
//	double temp;
//	int ipair;
//	int table, index;
//	double radius_cut;

//	for (ii = 0; ii < inum; ii++) {
//		i = ilist[ii];
//		itype = type[i];
//		itag = tag[i];
//		ix = x[i][0];
//		iy = x[i][1];
//		iz = x[i][2];
//		ivx = v[i][0];
//		ivy = v[i][1];
//		ivz = v[i][2];
//		fijn = 0.0;
//		fijt = 0.0;
//		jlist = firstneigh[i];
//		jnum = numneigh[i];
//		// loop for all neighbor particles
//		for (jj = 0; jj < jnum; jj++) {
//			j = jlist[jj];
//			jtype = type[j];
//			jtag = tag[j];
//			if (setflag[itype][jtype]) {
//				tags2str(str, itag, jtag);
// 				ipair = pair_list->find_hash(str);
//				if (ipair == -1) {
//					error->all(FLERR, "The pair does not exsit");
//				}
//				table = TABLE(ipair);
//				index = INDEX(ipair);
//				// rij = ri - rj (relative position)
//				rijx = ix - x[j][0];
//				rijy = iy - x[j][1];
//				rijz = iz - x[j][2];
//				rijsq = (rijx*rijx + rijy*rijy + rijz*rijz);
//				// withint the cutoff range
//				if (rijsq < cutsq[itype][jtype]) {
//					rij = sqrt(rijsq);
//					if (rij < EPSILON) continue;    // rij cannot be 0.0 in DEM systems
//					
//					if (radius_flag == 1) {
//						radius_cut = radius[i] + radius[j];
//						if (rij >= radius_cut) {
//							hash_pair[table][index].ivalues[0] = 0;
//							hash_pair[table][index].dvalues[0] = 0.0;
//							hash_pair[table][index].dvalues[1] = 0.0;
//							hash_pair[table][index].dvalues[2] = 0.0;
//							continue;
//						}
//					}
//
//					rij_inv = 1.0/rij;
//					// nij = rij/|rij|
//					nx = rijx * rij_inv;              
//					ny = rijy * rij_inv;
//					nz = rijz * rij_inv;
//					n[0] = nx;
//					n[1] = ny;
//					n[2] = nz;
//
//					// overlap between two particles = rcut - (rij . nij) = rcut - |rij|
//					// may need to be changed in the future by replaicng with particles's radius
//					if (radius_flag) drijn = radius_cut - rij;
//				    else drijn = cut[itype][jtype] - rij;
//					drijnx = drijn * nx;
//					drijny = drijn * ny;
//					drijnz = drijn * nz;
//	
//					vijx = ivx - v[j][0];
//					vijy = ivy - v[j][1];
//					vijz = ivz - v[j][2];
//
//					if (radius_flag) {
//						// distance from the center of particle to the contact point
//						Li = (rij*rij + radius[i] * radius[i] - radius[j] * radius[j]) / (2 * rij);
//						Lj = rij - Li;
//						// vij = vi - vj + (wi cross Ri) - (wj cross Rj) (relative velocity)  
//						Vec_Cross_Prod_3D(omegainij, omega[i], n);  // Ri = -n  (Ri: from center of particle i to the contact point
//						Vec_Cross_Prod_3D(omegajnij, omega[j], n);  // but n is the vector from particle j to particle i
//
//						vijx += -Li*omegainij[0] - Lj*omegajnij[0];
//						vijy += -Li*omegainij[1] - Lj*omegajnij[1];
//						vijz += -Li*omegainij[2] - Lj*omegajnij[2];
//					}
//
//					// |vijn| = vij . nij
//					vijn = vijx * nx + vijy * ny + vijz * nz;
//				    // vijn = |vijn| . nij
//					vijnx = vijn*nx;
//					vijny = vijn*ny;
//					vijnz = vijn*nz;
//
//					// vijt = vij - (vij . nij) . nij
//					vijtx = vijx - vijnx;
//					vijty = vijy - vijny;
//					vijtz = vijz - vijnz;
//					vijt = sqrt(vijtx*vijtx + vijty*vijty + vijtz*vijtz);
//					if (fabs(vijt) < EPSILON) vijt_inv = 0.0;
//					else vijt_inv = 1.0/vijt;
//
//					// tij = tij/|tij|
//					tx = vijtx * vijt_inv;
//					ty = vijty * vijt_inv;
//					tz = vijtz * vijt_inv;
//
//					// calculate drijt
//					// first time to contact
//					contact_flag = hash_pair[table][index].ivalues[0];
//					if (contact_flag == 0) {
//						if (fabs(vijn) < EPSILON) temp = dt;
//						else temp = MIN(fabs(drijn/vijn), dt);
//						drijtx = vijtx * temp;
//						drijty = vijty * temp;
//						drijtz = vijtz * temp;
//						hash_pair[table][index].ivalues[0] = 1;
//					}
//					// during the same contact
//					else {
//						drijtx = hash_pair[table][index].dvalues[0];
//						drijty = hash_pair[table][index].dvalues[1];
//						drijtz = hash_pair[table][index].dvalues[2];
//						
//						// update the tang. disp. for the next time step
//						drijtx = drijtx + vijtx*dt;
//						drijty = drijty + vijty*dt;
//						drijtz = drijtz + vijtz*dt;
//
//						temp = drijtx*nx + drijty*ny + drijtz*nz;
//
//						drijtx = drijtx - temp*nx;
//						drijty = drijty - temp*ny;
//						drijtz = drijtz - temp*nz;
//					}
//
//					hash_pair[table][index].dvalues[0] = drijtx;
//					hash_pair[table][index].dvalues[1] = drijty;
//					hash_pair[table][index].dvalues[2] = drijtz;
//					
//					// fijn = -kn*delta_rijn - Cn*vijn
//					fijnx = kn[itype][jtype]*drijnx - Cn[itype][jtype]*vijnx;
//					fijny = kn[itype][jtype]*drijny - Cn[itype][jtype]*vijny;
//					fijnz = kn[itype][jtype]*drijnz - Cn[itype][jtype]*vijnz;
//
//					fijn = sqrt(fijnx*fijnx + fijny*fijny + fijnz*fijnz);
//					
//					// tangential force
//					// fijt = -kt*delta_rijt - Ct*vijt
//					fijtx = -kt[itype][jtype] * drijtx - Ct[itype][jtype] * vijtx;
//					fijty = -kt[itype][jtype] * drijty - Ct[itype][jtype] * vijty;
//					fijtz = -kt[itype][jtype] * drijtz - Ct[itype][jtype] * vijtz;
//
//					fijt = sqrt(fijtx*fijtx + fijty*fijty + fijtz*fijtz); 
//					temp = fabs(mu[itype][jtype]*fijn);
//					// |fijt| > mu*|fnij|
//					if (fijt > temp) {
//						fijtx = -temp * tx;
//						fijty = -temp * ty;
//						fijtz = -temp * tz;
//					}
//
//					f[i][0] += (fijnx + fijtx);
//					f[i][1] += (fijny + fijty);
//					f[i][2] += (fijnz + fijtz);
//
//					f[j][0] -= (fijnx + fijtx);
//					f[j][1] -= (fijny + fijty);
//					f[j][2] -= (fijnz + fijtz);
//
//					if (particle->torque_flag) {
//						fijt_a[0] = fijtx;
//						fijt_a[1] = fijty;
//						fijt_a[2] = fijtz;
//						Vec_Cross_Prod_3D(torqueij, n, fijt_a);
//
//						torque[i][0] += -Li * torqueij[0];
//						torque[i][1] += -Li * torqueij[1];
//						torque[i][2] += -Li * torqueij[2];
//
//						torque[j][0] += -Lj * torqueij[0];
//						torque[j][1] += -Lj * torqueij[1];
//						torque[j][2] += -Lj * torqueij[2];
//					}
//					
//				} // if (rij < cut[itype][jtype])
//				else {
//					hash_pair[table][index].ivalues[0] = 0;
//					hash_pair[table][index].dvalues[0] = 0.0;
//					hash_pair[table][index].dvalues[1] = 0.0;
//					hash_pair[table][index].dvalues[2] = 0.0;
//				} // if (rij >= cut[itype][jtype])
//			}
//		} // for (jj = 0; jj < nlocal; jj++)
//	} // for (i = 0; i < nlocal; i++)
	
	cudaError_t error_t;
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	gpuComputedemforce << < int(nlocal + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >
	 (particle->devCoordX, particle->devCoordY, particle->devCoordZ,
	  particle->devVeloX, particle->devVeloY, particle->devVeloZ,
	  neighbor->devPairtable, neighbor->devNumneigh, 
//	  neighbor->devContactFlag, neighbor->devDrijtx, neighbor->devDrijty, neighbor->devDrijtz,
	  particle->devMass, particle->devRadius, particle->devType,
	  particle->devMask, nlocal, dt, radius_flag, particle->torque_flag, devSetflag, 
	  devCutsq, devKn, devKt, devCn, devCt, devMu,
//	  particle->devOmegaX, particle->devOmegaY, particle->devOmegaZ,
//	  particle->devTorqueX, particle->devTorqueY, particle->devTorqueZ,
	  particle->devForceX, particle->devForceY, particle->devForceZ);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairDEMLsd::set_pairlist()
{
	int ibucket;
	int ipair;
	int table, index;
	char str[BYTES];

	int i, j, ii, jj, inum, jnum, itag, jtag, itype, jtype;
	int *ilist;                                    
	int *jlist;
	int **list;
	int **firstneigh;
	int *numneigh;

	ilist = neighbor->neighlist->ilist;
	inum = neighbor->neighlist->inum;
	firstneigh = neighbor->neighlist->firstneigh;
	numneigh = neighbor->neighlist->numneigh;

	int *tag = particle->tag;
	int *type = particle->type;
	
	// prepare for rebuilding hash_pair
	pair_list->set_hash();

	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		itype = type[i];
		itag = tag[i];
		jlist = firstneigh[i];
		jnum = numneigh[i];
		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];
			jtype = type[j];
			jtag = tag[j];
			if (setflag[itype][jtype] == 1) {
				tags2str(str, itag, jtag);
				ipair = pair_list->find_hash(str);
				if (pair_list->nbuilds_total == 1 && ipair == -1) ipair = pair_list->insert_hash(str);
				else ipair = pair_list->set_pair(str);
			}
		}
	}

	// only get pair info. from neighbor processors
	// when nprocs > 1 and not the first build
	if (pair_list->nbuilds_total > 1 && parallel->nprocs > 1) pair_list->exchange_pair();
}

/* ----------------------------------------------------------------------
   str = itag-jtag
------------------------------------------------------------------------- */

void PairDEMLsd::tags2str(char *str, int itag, int jtag) 
{
	sprintf(str, "%d-%d", itag, jtag);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairDEMLsd::init_one(int i, int j)
{
	force->type2pair[i][j] = pair_id;

	cut[j][i] = cut[i][j];
	cutsq[j][i] = cutsq[i][j];
	kn[j][i] = kn[i][j];
	Cn[j][i] = Cn[i][j];
	kt[j][i] = kt[i][j];
	Ct[j][i] = Ct[i][j];
	mu[j][i] = mu[i][j];
	setflag[j][i] = setflag[i][j];
	force->type2pair[j][i] = force->type2pair[i][j];
}

/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairDEMLsd::set_style(int narg, char **arg)
{
	if (narg != 3) error->all(FLERR,"Illegal PairDEMLsd command");

	int n = strlen(arg[0]) + 1;
	style = new char[n];
	strcpy(style,arg[0]);

	e = atof(arg[1]);           // restitution coefficient

	cut_global = atof(arg[2]);

	// reset cutoffs that have been explicitly set

	if (allocated) {
		int i, j;
		for (i=1; i <= particle->ntypes; i++) 
		for (j = i+1; j <= particle->ntypes; j++) {
			if (setflag[i][j]) cut[i][j] = cut_global;
		}
	}
}

/* ----------------------------------------------------------------------
                Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairDEMLsd::set_coeff(int narg, char **arg)
{
	int ilo,ihi,jlo,jhi;

	if (narg != 8) error->all(FLERR, "Incorrect args for pair_coeff of pair_style dem/lsd"); 
	if(!allocated) allocate();

	// find the lower and upper bound of type set in the pair_coeff
    force->bounds(arg[0],particle->ntypes,ilo,ihi);
    force->bounds(arg[1],particle->ntypes,jlo,jhi);

	double cut_one = cut_global;
	cut_one = atof(arg[7]);
	
	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo,i); j <= jhi; j++) {
			kn[i][j] = atof(arg[2]);
			Cn[i][j] = atof(arg[3]);
			kt[i][j] = atof(arg[4]);
			Ct[i][j] = atof(arg[5]);
			mu[i][j] = atof(arg[6]);
			cut[i][j] = cut_one;
			cutsq[i][j] = cut_one * cut_one;
			//if (cut_max < cut_one) cut_max = cut_one;
			setflag[i][j] = 1;        // "1" means this pair has been setup 
			hostSetflag[i * TYPEMAX + j] = setflag[i][j];
			hostCutsq[i * TYPEMAX + j] = cutsq[i][j];
			hostKn[i * TYPEMAX + j] = kn[i][j];
			hostKt[i * TYPEMAX + j] = kt[i][j];
			hostCn[i * TYPEMAX + j] = Cn[i][j];
			hostCt[i * TYPEMAX + j] = Ct[i][j];
			hostMu[i * TYPEMAX + j] = mu[i][j];
			hostSetflag[j * TYPEMAX + i] = setflag[i][j];
			hostCutsq[j * TYPEMAX + i] = cutsq[i][j];
			hostKn[j * TYPEMAX + i] = kn[i][j];
			hostKt[j * TYPEMAX + i] = kt[i][j];
			hostCn[j * TYPEMAX + i] = Cn[i][j];
			hostCt[j * TYPEMAX + i] = Ct[i][j];
			hostMu[j * TYPEMAX + i] = mu[i][j];
			count++;
		}
	}
	if (count == 0) {
		error->all(FLERR,"Incorrest args for pair coefficients");
	}
	
    // setup for GPU parameters
    cudaError_t cudaStatus;
   
    cudaMalloc(&devSetflag, TYPEMAX * TYPEMAX * sizeof(int));
    cudaMalloc(&devCutsq, TYPEMAX * TYPEMAX * sizeof(double));
	cudaMalloc(&devKn, TYPEMAX * TYPEMAX * sizeof(double));
	cudaMalloc(&devKt, TYPEMAX * TYPEMAX * sizeof(double));
	cudaMalloc(&devCn, TYPEMAX * TYPEMAX * sizeof(double));
	cudaMalloc(&devCt, TYPEMAX * TYPEMAX * sizeof(double));
	cudaMalloc(&devMu, TYPEMAX * TYPEMAX * sizeof(double));
	   
    cudaMemcpy(devSetflag, hostSetflag, TYPEMAX * TYPEMAX * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(devKn, hostKn, TYPEMAX * TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devKt, hostKt, TYPEMAX * TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCn, hostCn, TYPEMAX * TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCt, hostCt, TYPEMAX * TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devMu, hostMu, TYPEMAX * TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
    cudaStatus = cudaMemcpy(devCutsq, hostCutsq, TYPEMAX * TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
  
}

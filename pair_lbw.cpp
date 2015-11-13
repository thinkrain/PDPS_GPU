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

#include "error.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair_lbw.h"
#include "parallel.h"
#include "particle.h"
#include "phy_const.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;
using namespace PhyConst;

#define DELTA 1
#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairLBW::PairLBW(PDPS *ps) : Pair(ps)
{
	phi = 0.0;

	allocated = 0;
	drag_flag = 0;

	newton_pair = 1;
}

/* ---------------------------------------------------------------------- */

PairLBW::~PairLBW()
{
	if (allocated) {
		memory->destroy(cut);
		memory->destroy(cutsq);
		memory->destroy(setflag);
	}
}

/* ---------------------------------------------------------------------- */

void PairLBW::allocate()
{
	allocated = 1;

	int n = particle->ntypes;

	cut = memory->create(cut, n + 1, n + 1, "PairLBW: cut");
	cutsq = memory->create(cutsq, n + 1, n + 1, "PairLBW: cutsq");
	setflag = memory->create(setflag, n + 1, n + 1, "PairLBW: setflag");

	// initialize setflag for pair i and j
	for (int i = 0; i <= n; i++)
	for (int j = 0; j <= n; j++) {
		setflag[i][j] = 0;
	}
}

/* ----------------------------------------------------------------------
Compute force for all paritcles
------------------------------------------------------------------------- */

void PairLBW::compute(int eflag, int vflag)
{
	int i, j, ii, jj, inum, jnum, itype, jtype;
	double ix, iy, iz, ivx, ivy, ivz;                    // pos and vel for particle i
	double rij, rijsq, rij_inv, rijx, rijy, rijz;        // relative position
	double drijn, drijnx, drijny, drijnz;                // normal displacement
	double drijtx, drijty, drijtz;                       // tangential displacement
	double vijx, vijy, vijz;                             // relative velocity: vij = vj - vi
	double vijn, vijnx, vijny, vijnz;                    // relative velocity along normal direction
	double vijtx, vijty, vijtz;                          // relative velocity along tangential direction
	double fvijnx, fvijny, fvijnz;                       // normal viscous force
	double fvijtx, fvijty, fvijtz;                       // tangential viscous force	
	double nx, ny, nz;                                   // unit vector along normal direction
	double tx, ty, tz;                                   // unit vector along tangential direction
	double fpair, evdwl;

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	int *type = particle->type;
	int radius_flag = particle->radius_flag;
	double *radius = particle->radius;
	int numlist;
	int nlocal = particle->nlocal;
	int nall = nlocal + particle->nghost;

	int *ilist;
	int *jlist;
	int **firstneigh;
	int *numneigh;

	double dt = update->dt;
	double dtinvsqrt = 1.0 / sqrt(dt);   // inverse of time step

	ilist = neighbor->neighlist->ilist;
	inum = neighbor->neighlist->inum;
	firstneigh = neighbor->neighlist->firstneigh;
	numneigh = neighbor->neighlist->numneigh;

	evdwl = 0.0;

	if (eflag || vflag) ev_setup(eflag, vflag);

	// loop for all local owned particles 
	double sij, Sij_half, R12, scut;
	double temp;
	double f1, f2, f3, f4;
	liquid_volume = liquid_volume_total = 0.0;
	ncollisions = ncollisions_total = 0;
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		itype = particle->type[i];
		ix = x[i][0];
		iy = x[i][1];
		iz = x[i][2];
		ivx = v[i][0];
		ivy = v[i][1];
		ivz = v[i][2];

		fvijnx = fvijny = fvijnz = 0.0;
		fvijtx = fvijty = fvijtz = 0.0;
		fpair = 0.0;
		jlist = firstneigh[i];
		jnum = numneigh[i];
		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];
			jtype = type[j];
			if (setflag[itype][jtype]) {
				rijx = ix - x[j][0];
				rijy = iy - x[j][1];
				rijz = iz - x[j][2];
				rijsq = (rijx*rijx + rijy*rijy + rijz*rijz);
				rij = sqrt(rijsq);
				sij = rij - radius[i] - radius[j];
				R12 = 2 * radius[i] * radius[j] / (radius[i] + radius[j]);
				scut = (1 + phi / 2.0) * R12 * Vpcb_cube_root;
				// withint the seperation cutoff range
				if (sij < scut) {
					if (sij < sij_min) sij = sij_min;    // rij can be 0.0 in DPD systems
					Sij_half = sij / 2 / R12;

					ncollisions += 1;
					liquid_volume += Vpcb * R12 * R12 * R12;  // local liquid volume

					f1 = ff1();
					f2 = ff2();
					f3 = ff3();
					f4 = ff4();
					temp = log(Sij_half / Vpcb_square_root);
					fpair = -2 * PI* R12 * gamma * exp(f1 - f2 * exp(f3*temp + f4*temp*temp));

					if (rij < EPSILON) rij_inv = 0.0;
					else rij_inv = 1.0 / rij;

					nx = rijx * rij_inv;
					ny = rijy * rij_inv;
					nz = rijz * rij_inv;

					if (drag_flag == 1) {
						vijx = ivx - v[j][0];
						vijy = ivy - v[j][1];
						vijz = ivz - v[j][2];

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

						fvijnx = -6 * PI * mu * vijnx * R12 * R12 / sij;
						fvijny = -6 * PI * mu * vijny * R12 * R12 / sij;
						fvijnz = -6 * PI * mu * vijnz * R12 * R12 / sij;
						
						fvijtx = -6 * PI * mu * R12 * vijtx * (8.0 / 15 * log(R12 / sij) + 0.9588);
						fvijty = -6 * PI * mu * R12 * vijty * (8.0 / 15 * log(R12 / sij) + 0.9588);
						fvijtz = -6 * PI * mu * R12 * vijtz * (8.0 / 15 * log(R12 / sij) + 0.9588);
					}

					f[i][0] += nx*fpair;
					f[i][1] += ny*fpair;
					f[i][2] += nz*fpair;

					// add pair force to the other particle

					f[j][0] -= nx*fpair;
					f[j][1] -= ny*fpair;
					f[j][2] -= nz*fpair;

					if (drag_flag) {
						f[i][0] += fvijnx + fvijtx;
						f[i][1] += fvijny + fvijty;
						f[i][2] += fvijnz + fvijtz;

						f[j][0] -= (fvijnx + fvijtx);
						f[j][1] -= (fvijny + fvijty);
						f[j][2] -= (fvijnz + fvijtz);
					}
				}
			}
		} // for (jj = 0; jj < nlocal; jj++)
	} // for (i = 0; i < nlocal; i++)

	MPI_Allreduce(&ncollisions, &ncollisions_total, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&liquid_volume, &liquid_volume_total, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (procid == 0 && file) {
		if (update->ntimestep % nevery == 0) {
			fprintf(file, "%d %d %g\n", update->ntimestep, ncollisions, liquid_volume_total);
			fflush(file);
		}
	}
}

/* ----------------------------------------------------------------------
init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairLBW::init_one(int i, int j)
{
	force->type2pair[i][j] = pair_id;

	cut[j][i] = cut[i][j];
	cutsq[j][i] = cutsq[i][j];
 
	setflag[j][i] = setflag[i][j];
	force->type2pair[j][i] = force->type2pair[i][j];
}

/* ----------------------------------------------------------------------
Setting for pair_style command
------------------------------------------------------------------------- */

void PairLBW::set_style(int narg, char **arg)
{
	if (narg < 6) error->all(FLERR, "Illegal PairLBW command");

	int n = strlen(arg[0]) + 1;
	style = new char[n];
	strcpy(style, arg[0]);               // store pair style's name

	double R1, R2, R12;
	double Srup;

	gamma = atof(arg[1]);
	Vpcb = atof(arg[2]);
	ln_Vpcb = log(Vpcb);
	Vpcb_square_root = sqrt(Vpcb);
	Vpcb_cube_root = pow(Vpcb, 1.0 / 3.0);
	sij_min = atof(arg[3]);
	R1 = atof(arg[4]);
	R2 = atof(arg[5]);

	R12 = 2*R1*R2 / (R1 + R2);

	Srup = 0.5 * (1 + phi / 2.0) * R12 * Vpcb_cube_root;

	cut_global = R1 + R2 + 2*Srup;       // store the gloabl cutoff

	int iarg = 6;
	while (iarg < narg) {
		if (!strcmp(arg[iarg], "file")) {
			n = strlen(arg[iarg+1]) + 1;
			fname = new char[n];
			strcpy(fname, arg[iarg+1]);
			if (procid == 0) {
				file = fopen(fname, "w");
				if (file == NULL) {
					char str[128];
					sprintf(str, "Cannot open analyze ave/space file %s", arg[iarg + 1]);
					error->one(FLERR, str);
				}
				fprintf(file, "step ncollisions liquid_volume\n");
			}
			iarg += 2;
		}
		else if (!strcmp(arg[iarg], "every")) {
			nevery = atoi(arg[iarg+1]);
			iarg += 2;
		}
		else if (!strcmp(arg[iarg], "drag")) {
			drag_flag = 1;
			mu = atof(arg[iarg+1]);
			iarg += 2;
		}
		else error->all(FLERR, "Illegal pair lb/willet option");
	}

	// reset cutoffs that have been explicitly set

	if (allocated) {
		int i, j;
		for (i = 1; i <= particle->ntypes; i++)
		for (j = i + 1; j <= particle->ntypes; j++) {
			if (setflag[i][j]) cut[i][j] = cut_global;
		}
	}
}

/* ----------------------------------------------------------------------
   Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairLBW::set_coeff(int narg, char **arg)
{
	int ilo, ihi, jlo, jhi;

	if (narg != 4) error->all(FLERR, "Incorrect args for pair_coeff of pair_style dpd");
	if (!allocated) allocate();

	// find the lower and upper bound of type set in the pair_coeff
	force->bounds(arg[0], particle->ntypes, ilo, ihi);
	force->bounds(arg[1], particle->ntypes, jlo, jhi);

	double R1, R2, R12;
	double Srup; 

	R1 = atof(arg[2]);
	R2 = atof(arg[3]);
	
	R12 = 2*R1*R2 / (R1 + R2);

	Srup = 0.5 * (1 + phi / 2.0) * R12 * Vpcb_cube_root;

	double cut_one;
	cut_one = R1 + R2 + 2 * Srup;

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			cut[i][j] = cut_one;
			cutsq[i][j] = cut_one * cut_one;
			setflag[i][j] = 1;        // "1" means this pair has been setup
			count++;
		}
	}

	if (count == 0) {
		error->all(FLERR, "Incorrest args for pair coefficients");
	}
}

/* ----------------------------------------------------------------------
   Check the paper for details:
   "2010 Mixing characteristics of wet granular matter in a bladed mixer"
   "2000 Capillary bridges between two sphere bodies"
   The following function is only valid when phi < 50 degree and V_pcb < 0.1
------------------------------------------------------------------------- */

double PairLBW::ff1()
{
	double res;

	res = (-0.44507 + 0.050832*phi - 1.1466*phi*phi)
		+ (-0.1119 - 0.000411*phi - 0.1490*phi*phi)*ln_Vpcb
		+ (-0.012101 - 0.0036456*phi - 0.01255*phi*phi)*ln_Vpcb*ln_Vpcb
		+ (-0.0005 - 0.0003505*phi - 0.00029076*phi*phi)*ln_Vpcb*ln_Vpcb*ln_Vpcb;

	return res;
}

/* ---------------------------------------------------------------------- */

double PairLBW::ff2()
{
	double res;

	res = (1.9222 - 0.57473*phi - 1.2918*phi*phi)
		+ (-0.0668 - 0.1201*phi - 0.22574*phi*phi)*ln_Vpcb
		+ (-0.0013375 - 0.0068988*phi - 0.01137*phi*phi)*ln_Vpcb*ln_Vpcb;

	return res;
}

/* ---------------------------------------------------------------------- */

double PairLBW::ff3()
{
	double res;

	res = (1.268 - 0.01396*phi - 0.23566*phi*phi)
		+ (0.198 + 0.092*phi - 0.06418*phi*phi)*ln_Vpcb
		+ (0.02232 + 0.02238*phi - 0.009853*phi*phi)*ln_Vpcb*ln_Vpcb
		+ (0.0008585 + 0.001318*phi - 0.00053*phi*phi)*ln_Vpcb*ln_Vpcb*ln_Vpcb;

	return res;
}

/* ---------------------------------------------------------------------- */

double PairLBW::ff4()
{
	double res;

	res = (-0.010703 + 0.073776*phi - 0.34742*phi*phi)
		+ (0.03345 + 0.04543*phi - 0.09056*phi*phi)*ln_Vpcb
		+ (0.0018574 + 0.004456*phi - 0.006257*phi*phi)*ln_Vpcb*ln_Vpcb;

	return res;
}

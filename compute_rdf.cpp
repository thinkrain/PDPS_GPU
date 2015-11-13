/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "compute_rdf.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair.h"
#include "particle.h"
#include "phy_const.h"
#include "update.h"

using namespace PDPS_NS;
using namespace PhyConst;

/* ---------------------------------------------------------------------- */

ComputeRdf::ComputeRdf(PDPS *ps, int narg, char **arg) : Compute(ps, narg, arg)
{
	if (narg < 4 || (narg - 4) % 2) error->all(FLERR,"Illegal compute rdf command");
	
	array_flag = 1;

	ilo = ihi = jlo = jhi = NULL;
	type2field = NULL;
	array = NULL;
	typecount = NULL;
	icount = jcount = NULL;
	type2field = NULL;
	gr = NULL;
	grAll = NULL;

	int ntypes = particle->ntypes;

	nbins = atoi(arg[3]);                 
	if (nbins < 1) error->all(FLERR,"# of bins has to be an integer larger than 0");
	if (narg == 4) nfields = 1;
	else nfields = (narg - 4) / 2;

	ilo = new int[nfields];
	ihi = new int[nfields];
	jlo = new int[nfields];
	jhi = new int[nfields];

	// rdf flag for i,j pair   default = -1
	memory->create(type2field, ntypes+1, ntypes+1, "ComputeRdf: type2field");
	for (int i = 0; i <= ntypes; i++) {
		for (int j = 0; j <= ntypes; j++) {
			type2field[i][j] = -1;
		}
	}
	// if no pair is specified
	if (narg == 4) {
		ilo[0] = 1; ihi[0] = ntypes;
		jlo[0] = 1; jhi[0] = ntypes;
	}
	else {
		int ifield;
		ifield = 0;
		int iarg =4;         // start at 4th command
		while (iarg < narg) {
			// find the lower and upper bound of type set in the pair_coeff
			force->bounds(arg[iarg], particle->ntypes, ilo[ifield], ihi[ifield]);
			force->bounds(arg[iarg+1], particle->ntypes, jlo[ifield], jhi[ifield]);
			if (ilo[ifield] > ihi[ifield] || jlo[ifield] > jhi[ifield]) {
				error->all(FLERR, "Illegal pair index");
			}
			ifield++;
			iarg += 2;
		}
	}

	// flag = 1 if i,j pair is specified
	for (int ifield = 0; ifield < nfields; ifield++) {
		for (int i = ilo[ifield]; i <= ihi[ifield]; i++) 
		for (int j = jlo[ifield]; j <= jhi[ifield]; j++) {
			type2field[i][j] = ifield;
		}
	}

	size_array_rows = nbins;
	size_array_columns = 1 + 2*nfields;

	memory->create(array, nbins, 1+2*nfields, "ComputeRdf: array");
	typecount = new int[ntypes+1];
	icount = new int[nfields];
	jcount = new int[nfields];

	memory->create(gr, nfields, nbins, "ComputeRdf: gr");
	memory->create(grAll, nfields, nbins, "ComputeRdf: gr");
}

/* ---------------------------------------------------------------------- */

ComputeRdf::~ComputeRdf()
{
	delete[] ilo;
	delete[] ihi;
	delete[] jlo;
	delete[] jhi;

	ilo = ihi = jlo = jhi = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeRdf::init()
{
	int i, ifield;

	int *mask = particle->mask;
	int *type = particle->type;
	int ntypes = particle->ntypes;
	int nlocal = particle->nlocal;

	rcut_max = 0.0;
	int pair_id;
	for (int ifield = 0; ifield < nfields; ifield++) {
		for (int i = ilo[ifield]; i <= ihi[ifield]; i++) 
		for (int j = jlo[ifield]; j <= jhi[ifield]; j++) {
			pair_id = force->type2pair[i][j];
			rcut_max = MAX(rcut_max, force->pair[pair_id]->cut[i][j]);
		}
	}

	rbin = rcut_max / nbins;
	rbinsinv = 1.0 / rbin;

	// set 1st column to bin coords
	for (int i = 0; i < nbins; i++) {
		array[i][0] = (i + 0.5) * rbin;
	}

	// count # of particles for each type within this group
	for (i = 0; i <= ntypes; i++) typecount[i] = 0;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) typecount[type[i]]++;
	}

	// count icount and jcount 
	// only need to do this once since all particles are within one group
	for (ifield = 0; ifield < nfields; ifield++) {
		icount[ifield] = 0;
		for (i = ilo[ifield]; i <= ihi[ifield]; i++) icount[ifield] += typecount[i];
		jcount[ifield] = 0;
		for (i = jlo[ifield]; i <= jhi[ifield]; i++) jcount[ifield] += typecount[i];
	}

	int *temp = new int[nfields];
	MPI_Allreduce(icount, temp, nfields, MPI_INT, MPI_SUM, mworld);
	for (i = 0; i < nfields; i++) icount[i] = temp[i];
	MPI_Allreduce(jcount, temp, nfields, MPI_INT, MPI_SUM, mworld);
	for (i = 0; i < nfields; i++) jcount[i] = temp[i];
	delete[] temp;
	temp = NULL;
}

/* ----------------------------------------------------------------------
   Compute rdf 
------------------------------------------------------------------------- */

void ComputeRdf::compute_array()
{
	int i, ii, j, jj, ifield;
	int *ilist, *jlist, *numneigh, **firstneigh;
	int itype, inum, jtype, jnum;
	double xtmp, ytmp, ztmp;
	double delx, dely, delz, rij;
	int ibin;
	int pair_id;

	invoked_array = update->ntimestep;
	
	ilist = neighbor->neighlist->ilist;
	inum = neighbor->neighlist->inum;
	firstneigh = neighbor->neighlist->firstneigh;
	numneigh = neighbor->neighlist->numneigh;

	double **x = particle->x;
	int *type = particle->type;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;

	class Pair **pair = force->pair;

	// zero gr[][]
	for (i = 0; i < nfields; i++) 
	for (j = 0; j < nbins; j++) {
		gr[i][j] = 0.0;
		grAll[i][j] = 0.0;
	}

	// loop for particle it owns
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		if (!(mask[i] & groupbit)) continue;
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		itype = type[i];
		jlist = firstneigh[i];
		jnum = numneigh[i];
		// loop for the neighbor list
		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];
			if (!(mask[j] & groupbit)) continue;
			jtype = type[j];
			ifield = type2field[itype][jtype];
			if (ifield < 0) continue;
			
			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rij = sqrt(delx*delx + dely*dely + delz*delz);
			ibin = static_cast<int> (rij*rbinsinv);
			if (ibin >= nbins) continue;
			
			gr[ifield][ibin] += 1.0;

			pair_id = force->type2pair[itype][jtype];
			if (pair[pair_id]->newton_pair == 1) {
				gr[ifield][ibin] += 1.0;
			}
		} // for (jj = 0; jj < jnum; jj++)
	} // for (ii = 0; ii < inum; ii++) 

	// sum g(r) across procs
	MPI_Allreduce(gr[0], grAll[0], nfields*nbins, MPI_DOUBLE, MPI_SUM, mworld);

	double constant, tmp, nums, rlower, rupper;
	double Vol_box, Vol_local;

	if (domain->dim = 3) {
		Vol_box = domain->xle*domain->yle*domain->zle;
		constant = 4.0 * PI / 3.0;  
		for (ifield = 0; ifield < nfields; ifield++) {
			nums = 0.0;
			for (ibin = 0; ibin < nbins; ibin++) {
				rlower = ibin * rbin;
				rupper = (ibin+1) * rbin;
				Vol_local = constant * (rupper*rupper*rupper - rlower*rlower*rlower);
				if (icount[ifield] == 0 || jcount[ifield] == 0) {
					gr[ifield][ibin] = 0.0;  // not quite sure if it is ncessary
					tmp = 0.0;
				}
				else {
					tmp = jcount[ifield] * Vol_local / Vol_box;
					gr[ifield][ibin] = grAll[ifield][ibin] / (icount[ifield]*tmp);
				}
				nums += gr[ifield][ibin] * tmp;
				array[ibin][1+2*ifield] = gr[ifield][ibin];
				array[ibin][2+2*ifield] = nums;
			}
		}
	} // if (domain->dim = 3)
	else {
		Vol_box = domain->xle*domain->yle;
		for (ifield = 0; ifield < nfields; ifield++) {
			nums = 0.0;
			for (ibin = 0; ibin < nbins; ibin++) {
				rlower = ibin * rbin;
				rupper = (ibin+1) * rbin;
				Vol_local = PI * (rupper*rupper - rlower*rlower);
				if (icount[ifield] == 0 || jcount[ifield] == 0) {
					gr[ifield][ibin] = 0.0;  // not quite sure if it is ncessary
				}
				else {
					tmp = jcount[ifield] * Vol_local / Vol_box;
					gr[ifield][ibin] = grAll[ifield][ibin] / (icount[ifield]*tmp);
				}
			}
			nums += gr[ifield][ibin] * tmp;
			array[ibin][1+2*ifield] = gr[ifield][ibin];
			array[ibin][2+2*ifield] = nums;
		}
	} // if (domain->dim == 2)
}

/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "error.h"
#include "lattice.h"
#include "memory.h"
#include "region.h"

using namespace PDPS_NS;

#define DELTA 4         // For reallocation of memory, if the size of the array need to be expanded 

enum{SPACING, CUSTOM};

/* ---------------------------------------------------------------------- */

Lattice::Lattice(PDPS *ps, int narg, char **arg) : Pointers(ps)
{
	//cell_coord = NULL;
	cell_num = NULL;

	if (strcmp(arg[1], "spacing") == 0) {
		rid = domain->find_region(arg[0]);
		if (rid == -1) error->all(FLERR, "Cannot find the region");
		region = domain->regions[rid];
		style = SPACING;
		rcell = 0.5*atof(arg[2]);
		setup_cell();
	}
	else error->all(FLERR, "Illegal lattice style");
}

/* ---------------------------------------------------------------------- */

Lattice::~Lattice()
{
	delete[] cell_num;
	cell_num = NULL;

	//memory->destroy(cell_coord);
}

/* ---------------------------------------------------------------------- */

void Lattice::setup_cell()
{
	int i, j, k;
	int dim = domain->dim;

	extent_lo = region->extent_lo;
	extent_hi = region->extent_hi;
	extent_le = region->extent_le;

	double coord;

	// create cells in the global domain
	nc[0] = static_cast <int> (extent_le[0]/rcell);
	nc[1] = static_cast <int> (extent_le[1]/rcell);
	if (dim == 3) nc[2] = static_cast <int> (extent_le[2]/rcell);
	else nc[2] = 1;

	// make sure there is one bin if rneigh > box size
	for (i = 0; i < 3; i++) {
		if (nc[i] == 0) nc[i] = 1;
	}

	ncxy = nc[0]*nc[1];
	ncxyz = ncxy*nc[2]; 

	cell_num = new int[ncxyz];
	for (i = 0; i < ncxyz; i++) {
		cell_num[i] = 0;
	}

	// update cell length, so that the allocation will be more uniform
	for (i = 0; i < 3; i++) cle[i] = extent_le[i] / nc[i];

	// allocate
	//memory->create(cell_coord, ncxyz, 3, "Lattice: cell_coord");

	/* int cid; 
	for (k = 0; k < nc[2]; k++)
	for (j = 0; j < nc[1]; j++)
	for (i = 0; i < nc[0]; i++) {
		cid = k*ncxy + j*nc[0] + i;
		cell_coord[cid][0] = i * cle[0] + extent_lo[0];
		cell_coord[cid][1] = j * cle[1] + extent_lo[1];
		cell_coord[cid][2] = k * cle[2] + extent_lo[2];
	}*/
}

/* ----------------------------------------------------------------------
   coordinate to cell id
------------------------------------------------------------------------- */

int Lattice::coord2cell(double *x, int &icx, int &icy, int &icz)
{
	int ic[3];

	for (int i = 0; i < 3; i++) {
		ic[i] = static_cast<int> ((x[i] - extent_lo[i])/cle[i]);
		ic[i] = MIN(ic[i],nc[i]-1);
	}

	if (domain->dim == 2) {
		ic[2] = 0;
	}

	icx = ic[0];
	icy = ic[1];
	icz = ic[2];

	int c = ic[2]*ncxy + ic[1]*nc[0] + ic[0];
	return c;
}

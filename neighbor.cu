/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */


#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "error.h"
#include "force.h"
#include "domain.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "parallel.h"
#include "particle.h"
#include "pair.h"
#include "update.h"
#include "group.h"
#include "timer.h"

#include "pdps_cuda.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "device_functions.h"

using namespace PDPS_NS;

#define RQDELTA 1
#define EXDELTA 1

#define coeff_neigh 0.5
#define coeff_delta 0.1

#define LB_FACTOR 1.5
#define SMALL 1.0e-6
#define BIG 1.0e20
#define CUT2BIN_RATIO 100


enum{SINGLE,MULTI};     // also in neigh_list.cpp

//#define NEIGH_LIST_DEBUG 1

/* ---------------------------------------------------------------------- */
__global__ void gpulinklist(int *devLinked_list, int *devHead, int *devSubclo, int *devSubnc,
	double *devBoxhi, double *devBoxlo, double *devCle, int *devNc, const int nall,
	const int subncxyz, const int subncxy, double *devCoordX, double *devCoordY, double *devCoordZ){


	// Reset the headers to -1
	for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < subncxyz; i+=blockDim.x * gridDim.x) 
		devHead[i] = -1;
	__syncthreads();

	__shared__ double boxhi[3];
	__shared__ double boxlo[3];
	__shared__ double cle[3];
	__shared__ int nc[3];
	__shared__ int subclo[3];
	__shared__ int subnc[3];
	for (int dim = 0; dim < 3; dim++){
		boxhi[dim] = devBoxhi[dim];
		boxlo[dim] = devBoxlo[dim];
		cle[dim] = devCle[dim];
		nc[dim] = devNc[dim];
		subclo[dim] = devSubclo[dim];
		subnc[dim] = devSubnc[dim];
	}
	__syncthreads();
	// Scan particles to construct headers, head, & linked lists, linked_list
	// Scan ghost particles first, so they will appear at the end of the neighbor list
	// Scan in an reverse order, so they will apear in a correct order in the neighbor list 

	for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < subncxyz; i += blockDim.x * gridDim.x){
		for (int pid = 0; pid <  nall; pid++){
			if (pid >= nall)
				break;
			int ic[3];
			for (int dim = 0; dim < 3; dim++) {
				if (devCoordX[pid] >= boxhi[dim])
					ic[dim] = static_cast<int> ((devCoordX[pid] - boxhi[dim]) / cle[dim]) + nc[dim];
				else if (devCoordX[pid] >= boxlo[dim]) {
					ic[dim] = static_cast<int> ((devCoordX[pid] - boxlo[dim]) / cle[dim]);
					ic[dim] = MIN(ic[dim], nc[dim] - 1);
				}
				else
					ic[dim] = static_cast<int> ((devCoordX[pid] - boxlo[dim]) / cle[dim]) - 1;
			}

			for (int dim = 0; dim < 3; dim++)
				ic[dim] -= subclo[dim];

			int c = ic[2] * subncxy + ic[1] * subnc[0] + ic[0];

			if (c == i){
				// Link to the previoius occupant (or -1 if you are the 1st)
				devLinked_list[pid] = devHead[c];
				// The last one goes to the header
				devHead[c] = pid;
			}

		}
	}
	


}
 
//	build the whole neighbor list, not half
__global__  void gpuneighbuild(int *devLinked_list, int *devHead, int *devNumneigh, int *devPairtable, int *devCoffsets, const int devNoffsets, 
	int *devSubclo, int *devSubnc, double *devBoxhi, double *devBoxlo, double *devCle, int *devNc, double *devRneighsq, int *devMask,
	const int nlocal, const int subncxyz, const int subncxy, double *devCoordX, double *devCoordY, double *devCoordZ, int *devType, const int slave_bit){

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int ic[3];
	int c1, j, itype, jtype, dim;
	double rijsq, dx, dy, dz;

	__shared__ double boxhi[3];
	__shared__ double boxlo[3];
	__shared__ double cle[3];
	__shared__ int nc[3];
	__shared__ int subclo[3];
	__shared__ int subnc[3];
	__shared__ double rneighsq[TYPEMAX * TYPEMAX];
	for (dim = 0; dim < 3; dim++){
		boxhi[dim] = devBoxhi[dim];
		boxlo[dim] = devBoxlo[dim];
		cle[dim] = devCle[dim];
		nc[dim] = devNc[dim];
		subclo[dim] = devSubclo[dim];
		subnc[dim] = devSubnc[dim];
	}
	for (itype = 0; itype < TYPEMAX; itype++){
		for (jtype = 0; jtype < TYPEMAX; jtype++)
			rneighsq[itype * TYPEMAX + jtype] = devRneighsq[itype * TYPEMAX + jtype];
	}
	__syncthreads();

	for (i = i; i < nlocal; i += blockDim.x * gridDim.x){
		if (devMask[i] & slave_bit)
			continue;
		int numneigh = 0;
		itype = devType[i];
		double xi = devCoordX[i];
		double yi = devCoordY[i];
		double zi = devCoordZ[i];
		dim = 0;
		//	X dimension
		if (xi >= boxhi[dim])
			ic[dim] = static_cast<int> ((xi - boxhi[dim]) / cle[dim]) + nc[dim];
		else if (xi >= boxlo[dim]) {
			ic[dim] = static_cast<int> ((xi - boxlo[dim]) / cle[dim]);
			ic[dim] = MIN(ic[dim], nc[dim] - 1);
		}
		else
			ic[dim] = static_cast<int> ((xi - boxlo[dim]) / cle[dim]) - 1;
		//  Y dimension
		dim = 1;
		if (yi >= boxhi[dim])
			ic[dim] = static_cast<int> ((yi - boxhi[dim]) / cle[dim]) + nc[dim];
		else if (yi >= boxlo[dim]) {
			ic[dim] = static_cast<int> ((yi - boxlo[dim]) / cle[dim]);
			ic[dim] = MIN(ic[dim], nc[dim] - 1);
		}
		else
			ic[dim] = static_cast<int> ((yi - boxlo[dim]) / cle[dim]) - 1;
		//  Z dimension
		dim = 2;
		if (zi >= boxhi[dim])
			ic[dim] = static_cast<int> ((zi - boxhi[dim]) / cle[dim]) + nc[dim];
		else if (zi >= boxlo[dim]) {
			ic[dim] = static_cast<int> ((zi - boxlo[dim]) / cle[dim]);
			ic[dim] = MIN(ic[dim], nc[dim] - 1);
		}
		else
			ic[dim] = static_cast<int> ((zi - boxlo[dim]) / cle[dim]) - 1;

		for (dim = 0; dim < 3; dim++)
			ic[dim] -= subclo[dim];

		int c = (ic[2] + 2) * subncxy + (ic[1] + 2) * subnc[0] + ic[0] + 2;

		//	search neighbor particles through cells, in both direction
		for (int direction = 0; direction < 2; direction++){
			for (int io = 0; io < devNoffsets; io++){
				if (direction == 0)
					c1 = c + devCoffsets[io];
				else
					c1 = c - devCoffsets[io];
				if (c == c1 && direction == 1) {
					continue;
				}
				else j = devHead[c1];

				while (j != -1){
					jtype = devType[j];
					// scan rest of particles in this cell c1
					if (c == c1) {
						// skip itself
						if (j == i){
							j = devLinked_list[j];
							continue;
						}
					}

					rijsq = 0.0;
					dx = xi - devCoordX[j];
					dy = yi - devCoordY[j];
					dz = zi - devCoordZ[j];
					rijsq += dx * dx;
					rijsq += dy * dy;
					rijsq += dz * dz;

					if (rijsq < rneighsq[itype * TYPEMAX + jtype]){
						devPairtable[i * NEIGHMAX + numneigh] = j;
						numneigh += 1;
					}

					j = devLinked_list[j];
				}	// j != -1

			}	// io < noffsets
		}	// direction < 2
		
		devNumneigh[i] = numneigh;


	}	// i < nlocal


}

//	build the whole neighbor list, not half
__global__  void gpuneighbuild2(int *devNumneigh, int *devPairtable, double *devRneighsq, const double Rcut, int *devMask,
	const int nlocal,double *devCoordX, double *devCoordY, double *devCoordZ, int *devType, const int slave_bit){

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j, itype, jtype;
	double rijsq, dx, dy, dz;

	__shared__ double rneighsq[TYPEMAX * TYPEMAX];

	for (itype = 0; itype < TYPEMAX; itype++){
		for (jtype = 0; jtype < TYPEMAX; jtype++)
			rneighsq[itype * TYPEMAX + jtype] = devRneighsq[itype * TYPEMAX + jtype];
	}
	__syncthreads();

	

	for (i = i; i < nlocal; i += blockDim.x * gridDim.x){
		if (devMask[i] & slave_bit)
			continue;
		int numneigh = 0;
		itype = devType[i];
		double xi = devCoordX[i];
		double yi = devCoordY[i];
		double zi = devCoordZ[i];
		
		for (int jj = i; jj < nlocal + i; jj++){
			if (jj >= nlocal)
				j = jj - nlocal;
			else
				j = jj;

			dx = xi - devCoordX[j];
			if (dx > -Rcut && dx < Rcut){
				dy = yi - devCoordY[j];
				if (dy > -Rcut && dy < Rcut){
					dz = zi - devCoordZ[j];
					if (dz > -Rcut && dz < Rcut){
						if (j == i)
							continue;
						rijsq = dx * dx + dy * dy + dz * dz;
						jtype = devType[j];
						if (rijsq < rneighsq[itype * TYPEMAX + jtype]){
							devPairtable[i * NEIGHMAX + numneigh] = j;
							numneigh++;
						}	//	rijsq
					}	//	dz
				}	//	dy
			}	// dx
		}	// j < nlocal
		devNumneigh[i] = numneigh;

	}	// i < nlocal


}

__global__  void gpucheck_distance(double *devCoordX, double *devCoordY, double *devCoordZ, double *devCoordXold, double *devCoordYold, double *devCoordZold,
							const int nlocal, const double delta, int *build_flag){
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	double dx, dy, dz, rsq;
	for (i = i; i < nlocal; i += blockDim.x * gridDim.x){
		dx = devCoordX[i] - devCoordXold[i];
		dy = devCoordY[i] - devCoordYold[i];
		dz = devCoordZ[i] - devCoordZold[i];
		rsq = dx * dx + dy * dy + dz * dz;
		if (rsq > delta)
			build_flag[0] = 1;
	}
}

// copy data between device
template<class TYPE> __global__ void gpuTransfer(
	TYPE* __restrict Out,
	TYPE* __restrict In,
	const int  n)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	for (i = i; i < n; i += blockDim.x * gridDim.x)
		Out[i] = In[i];
}

Neighbor::Neighbor(PDPS *ps) : Pointers(ps)
{
	procid = parallel->procid;

	head = NULL;
	head_slave = NULL;
	slave_linked = 0;
	linked_list = NULL;
    neighlist = NULL;
	rneigh = rneighsq = NULL;
	rneigh_type = rneighsq_type = NULL;
	x_old = NULL;

	ago = 0;
	boxcheck = 0;
	build_before = 0;
    delay = 1;
	dist_check = 1;
	every = 1;
	last_build = 0;
	nbuilds = nbuilds_total = 0;
	ndanger = ndanger_total = 0;
	neigh_flag = 0;
	nlocal_old = nmax_old = 0;
	rcut_max = 0.0;

	noffsets = 0;
	max_noffsets = 0;
	max_noffsets_multi = 0;
	coffsets = NULL;
	coffsets_multi = NULL;
	distsq_multi = NULL;
	noffsets_multi = NULL;

	subncxyz = subncxyz_max = 0;
	slave_flag = 0;

	// pgsize and neigh_threshold_one may be adjustable by user in the future
	pgsize = 100000;            // size of one page 
	onesize = 2000;             // max # of neighbor list allowed for one particle (default: 2000)

	neighlist = NULL;
	neighlist = new NeighList(ps,pgsize);

	devCoffsets = NULL;
	devLinked_list = NULL;
	devHead = NULL;
	devBoxle = NULL;
	devCle = NULL;
	devNc = NULL;

	devSubclo = NULL;
	devSubnc = NULL;
	devBoxhi = NULL;
	devBoxlo = NULL;

	devPairtable = NULL;
	devNumneigh = NULL;
}

/* ---------------------------------------------------------------------- */

Neighbor::~Neighbor()
{
	memory->destroy(head);
	memory->destroy(linked_list);
	memory->destroy(x_old);
	if (slave_flag)
		memory->destroy(head_slave);
	
	delete neighlist;
	neighlist = NULL;
	
	if (style == SINGLE) {
		memory->destroy(coffsets);
	}
	else if (style == MULTI) {
		if (coffsets_multi) {
			for (int i = 1; i <= particle->ntypes; i++) {
				memory->destroy(coffsets_multi[i]);
				memory->destroy(distsq_multi[i]);
			}
			delete[] coffsets_multi;
			coffsets_multi = NULL;
			delete[] distsq_multi;
			distsq_multi = NULL;
			delete[] noffsets_multi;
			noffsets_multi = NULL;
		}
	}

	memory->destroy(rneigh);
	memory->destroy(rneighsq);
	memory->destroy(rneigh_type);
	memory->destroy(rneighsq_type);

	cudaFree(devCoffsets);
	cudaFree(devLinked_list);
	cudaFree(devHead);
	cudaFree(devBoxle);
	cudaFree(devCle);
	cudaFree(devNc);
	cudaFree(devSubclo);
	cudaFree(devSubnc);
	cudaFree(devBoxhi);
	cudaFree(devBoxlo);
	cudaFree(devPairtable);
	cudaFree(devNumneigh);

}

/* ---------------------------------------------------------------------- */

void Neighbor::init()
{
	int i, j, k;
	nbuilds = ndanger = 0;

	if (delay > 0 && delay % every != 0) {
		error->all(FLERR, "delay should be multiple of every");
	}

	dim = domain->dim;
	nparticles = particle->nparticles;

	boxlo = domain->boxlo;
	boxhi = domain->boxhi;
	boxle = domain->boxle;
	sublo = domain->sublo;
	subhi = domain->subhi;
	suble = domain->suble;

	// sub lo/hi = bounding box of the local domain
	sublo_old[0] = sublo[0];
	sublo_old[1] = sublo[1];
	sublo_old[2] = sublo[2];
	subhi_old[0] = subhi[0];
	subhi_old[1] = subhi[1];
	subhi_old[2] = subhi[2];

	triggersq = 0.25*rskin*rskin;
	if (domain->box_change && (domain->xperiodic || domain->yperiodic ||
		(dim == 3 && domain->zperiodic))) {
		boxcheck = 1;
	}

	rcut_max = force->cut_max_global;
	rcut_min = force->cut_min_global;
	rneigh_max = rcut_max + rskin;
	rneigh_min = rcut_min + rskin;

	int n = particle->ntypes;
	if (rneigh == NULL) {
		memory->create(rneigh, n+1, n+1, "neighbor: rneigh");
		memory->create(rneighsq, n+1, n+1, "neighbor: rneighsq");
		memory->create(rneigh_type, n+1, "neighbor: rneigh_type");
		memory->create(rneighsq_type, n+1, "neighbor: rneighsq_type");
	}

	for (i = 0; i < n + 1; i++) {
		rneigh_type[i] = 0.0;
		rneighsq_type[i] = 0.0;
		for (j = 0; j < n + 1; j++) {
			rneigh[i][j] = 0.0;
			rneighsq[i][j] = 0.0;
		}
	}
	// Initialization all rneigh type related
	double rcut, rneigh_temp;
	int ipair;
	Pair **pair = force->pair;
	int npairs = force->npairs;
	for (i = 1; i <= n; i++) {
		rneigh_type[i] = 0.0;
		rneighsq_type[i] = 0.0;
		for (j = 1; j <= n; j++) {
			for (ipair = 0; ipair < npairs; ipair++) {
				if (pair[ipair]->setflag[i][j]) {
					rcut = pair[ipair]->cut[i][j];
					rneigh[i][j] = MAX(rneigh[i][j], rcut + rskin);
					rneighsq[i][j] = rneigh[i][j] * rneigh[i][j];
					rneigh_type[i] = MAX(rneigh_type[i], rneigh[i][j]);
					rneighsq_type[i] = MAX(rneighsq_type[i], rneighsq[i][j]);
				}
			}
		}
	}
	rneighsq_max = rneigh_max * rneigh_max;
	rneighsq_min = rneigh_min * rneigh_min;

	// the first time to add page
	if (rneigh_max != 0 && neighlist->maxpage == 0) {	
		neighlist->add_pages(1);
		neighlist->list_exist = 1;
	} 

	// choose neighbor build method
	if (style == SINGLE) {
		pair_build = &Neighbor::half_linkedlist;
		offsets_create = &Neighbor::offsets_half_single;
		
	}
	else if (style == MULTI) {
		pair_build = &Neighbor::half_multi;
		offsets_create = &Neighbor::offsets_half_multi;
	}

	double *hostRneighsq;
	hostRneighsq = (double *)malloc(TYPEMAX * TYPEMAX * sizeof(double));
	for (i = 1; i < n; i++){
		for (j = 1; j < n; j++)
			hostRneighsq[i * TYPEMAX + j] = rneighsq[i][j];
	}
	cudaMalloc(&devRneighsq, TYPEMAX * TYPEMAX * sizeof(double));
	cudaMemcpy(devRneighsq, hostRneighsq, TYPEMAX * TYPEMAX * sizeof(double), cudaMemcpyHostToDevice);
	int n_master;
	if (slave_flag)
		n_master = particle->nlocal - slave_num;
	else
		n_master = particle->nlocal;
	gpuTransfer << < int(n_master + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> > (particle->devCoordXold, particle->devCoordX, n_master);
	gpuTransfer << < int(n_master + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> > (particle->devCoordYold, particle->devCoordY, n_master);
	gpuTransfer << < int(n_master + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> > (particle->devCoordZold, particle->devCoordZ, n_master);

}

/* ---------------------------------------------------------------------- 
                            Settings
---------------------------------------------------------------------- */

void Neighbor::settings(int narg, char **arg)
{
	if (narg != 2) error->all(FLERR, "Illegal neighbor command");

	rskin = atof(arg[0]);         // store skin distance
	if (rskin <= 0) error->all(FLERR,"Skin value should be positive");

	if (strcmp(arg[1],"single") == 0) style = SINGLE;
	else if (strcmp(arg[1],"multi") == 0) style = MULTI;
	else error->all(FLERR, "Illegal neighbor command");
}

/* ---------------------------------------------------------------------- 
                        Set coefficients
---------------------------------------------------------------------- */

void Neighbor::modify_params(int narg, char **arg)
{
	int n = narg / 2;
	if (narg % 2 != 0) error->all(FLERR,"Illegal neigh_modify command");

	int iarg = 0;
	while (iarg < narg) {
		if (!strcmp(arg[iarg],"every")) every = atoi(arg[iarg+1]);
		else if (!strcmp(arg[iarg],"delay")) delay = atoi(arg[iarg+1]);
		else if (!strcmp(arg[iarg],"check")) {
			if (!strcmp(arg[iarg+1],"yes")) {
				dist_check = 1;
			}
			else {
				dist_check = 0;
			}
		}
		else if (!strcmp(arg[iarg],"one")) {
			onesize = atoi(arg[iarg+1]);
		}
		else error->all(FLERR, "Illegal neigh_modify command");
		iarg += 2;
	}

	if (every > delay || delay % every != 0) {
		error->all(FLERR,"delay should be multiple of every");
	}
}

/* ---------------------------------------------------------------------- 
   Setup cell for each extended sub-domain:
   1. cell must be setup in the global domain, so that the index of each
   particle and its image can be matched. For example, the domain is divided
   into 5 by 5 by 5. If a particle has index (1,0,2), one of its image (>boxxhi)
   should be exactly (6,0,2). Otherwise, during the build of neighbors, the pair
   information will not be consistent. 
---------------------------------------------------------------------- */

void Neighbor::setup_cells()
{
	int nmax = particle->nmax;
	allocate(nmax);

	int c;

	// ghost_sub lo/hi = bounding box of the local domain extended by parallel->rcghost
	
	double *rcghost = parallel->rcghost;
	// extend subdomain by rcghost
	ghost_sublo[0] = sublo[0] - rcghost[0];
	ghost_sublo[1] = sublo[1] - rcghost[1];
	ghost_sublo[2] = sublo[2] - rcghost[2];
	ghost_subhi[0] = subhi[0] + rcghost[0];
	ghost_subhi[1] = subhi[1] + rcghost[1];
	ghost_subhi[2] = subhi[2] + rcghost[2];

	double rcell;
	
	if (style == SINGLE) {
		rcell = 0.5 * rcut_max;         // default cell length is half of the maximum pair cut-off
	}
	else if (style == MULTI) {
		rcell = 0.5 * rcut_min;
	}
	
	// create cells in the global domain
	nc[0] = static_cast <int> (boxle[0]/rcell);
	nc[1] = static_cast <int> (boxle[1]/rcell);
	if (dim == 3) nc[2] = static_cast <int> (boxle[2]/rcell);
	else nc[2] = 1;
  
	// make sure there is one bin if rneigh > box size
	for (int i = 0; i < 3; i++) {
		if (nc[i] == 0) nc[i] = 1;
	}

	ncxy = nc[0]*nc[1];
	ncxyz = ncxy*nc[2]; 

	// update cell length, so that the allocation will be more uniform
	for (int i = 0; i < 3; i++) {
		cle[i] = boxle[i]/nc[i];
	}

	// find the lowerest and uppermost global index in the extended subdomain
	int subchi[3];
	double coord;

	for (int i = 0; i < 3; i++) {
		coord = ghost_sublo[i] - SMALL*boxle[i];
		subclo[i] = static_cast<int> ((coord - boxlo[i])/cle[i]);
		if (coord < boxlo[i]) subclo[i] = subclo[i] - 1;

		coord = ghost_subhi[i] + SMALL*boxle[i];
		subchi[i] = static_cast<int> ((coord - boxlo[i])/cle[i]);
	}
	if (dim == 2) {
		subclo[2] = 0;
		subchi[2] = 0;
	}

	// extend cell by 1 to insure cell offsets can cover all possible interactive cells
	// if 2d, only 1 cell in z

	subclo[0] -= 2;
	subchi[0] += 2;
	subclo[1] -= 2;
	subchi[1] += 2;

	if (dim == 3) {
		subclo[2] -= 2;
		subchi[2] += 2;
	}

	for (int i = 0; i < 3; i++) {
		subnc[i] = subchi[i] - subclo[i] + 1;
	}

	// total number of cells to be created in local domain
	subncxy = subnc[0] * subnc[1];
	subncxyz = subncxy * subnc[2];
	
	if (subncxyz > subncxyz_max) {
		subncxyz_max = subncxyz;
		memory->destroy(head);
		head = memory->create(head, subncxyz_max, "Neighbor: head");
		if (slave_flag)
			head_slave = memory->create(head_slave, subncxyz_max, "Neighbor: head_slave");
	}

	(this->*offsets_create)();
	//	Transfer the data to GPU
	cudaMalloc(&devCoffsets, noffsets * sizeof(int));
	cudaMalloc(&devBoxle, 3 * sizeof(double));
	cudaMalloc(&devCle, 3 * sizeof(double));
	cudaMalloc(&devNc, 3 * sizeof(int));
	cudaMemcpy(devCoffsets, coffsets, noffsets * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devBoxle, boxle, 3 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCle, cle, 3 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devNc, nc, 3 * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc(&devSubclo, 3 * sizeof(int));
	cudaMalloc(&devSubnc, 3 * sizeof(int));
	cudaMalloc(&devBoxhi, 3 * sizeof(double));
	cudaMalloc(&devBoxlo, 3 * sizeof(double));
	cudaMemcpy(devSubclo, subclo, 3 * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devSubnc, subnc, 3 * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devBoxhi, boxhi, 3 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devBoxlo, boxlo, 3 * sizeof(double), cudaMemcpyHostToDevice);

	cudaMalloc(&devLinked_list, (particle->nlocal + particle->nghost) * sizeof(int));
	cudaMalloc(&devHead, subncxyz * sizeof(int));
	cudaMalloc(&devNumneigh, particle->nlocal * sizeof(int));
	cudaMalloc(&devPairtable, particle->nlocal * NEIGHMAX * sizeof(int));

	hostPairtable = (int *)malloc(particle->nlocal * NEIGHMAX * sizeof(int));
	hostNumneigh = (int *)malloc(particle->nlocal * sizeof(int));
	hostCoordX = (double *)malloc((particle->nlocal + particle->nghost) * sizeof(double));
	hostCoordY = (double *)malloc((particle->nlocal + particle->nghost) * sizeof(double));
	hostCoordZ = (double *)malloc((particle->nlocal + particle->nghost) * sizeof(double));
	hostForceX = (double *)malloc((particle->nlocal + particle->nghost) * sizeof(double));
	hostForceY = (double *)malloc((particle->nlocal + particle->nghost) * sizeof(double));
	hostForceZ = (double *)malloc((particle->nlocal + particle->nghost) * sizeof(double));
	hostVeloX = (double *)malloc((particle->nlocal + particle->nghost) * sizeof(double));
	hostVeloY = (double *)malloc((particle->nlocal + particle->nghost) * sizeof(double));
	hostVeloZ = (double *)malloc((particle->nlocal + particle->nghost) * sizeof(double));

}

/* ----------------------------------------------------------------------
   compute closest distance between central bin (0,0,0) and bin (i,j,k)
------------------------------------------------------------------------- */

double Neighbor::cell_distance(int i, int j, int k)
{
	double delx,dely,delz;

	if (i > 0) delx = (i - 1)*cle[0];
	else if (i == 0) delx = 0.0;
	else delx = (i + 1)*cle[0];

	if (j > 0) dely = (j - 1)*cle[1];
	else if (j == 0) dely = 0.0;
	else dely = (j + 1)*cle[1];

	if (k > 0) delz = (k - 1)*cle[2];
	else if (k == 0) delz = 0.0;
	else delz = (k + 1)*cle[2];

	return (delx*delx + dely*dely + delz*delz);
}

/* ---------------------------------------------------------------------- 
   Build neighbor list
---------------------------------------------------------------------- */

void Neighbor::build()
{

	/*cudaError_t error_t;
	error_t = cudaMemcpy(hostPairtable, devPairtable, 1 * sizeof(int), cudaMemcpyDeviceToHost);
*/


	ago = 0;
	nbuilds++;
	nbuilds_total++;

	last_build = update->ntimestep;
	int nlocal = particle->nlocal;
	int nall = nlocal + particle->nghost;
	
	// box lo/hi = size of bbox of entire domain
	// store box dimension information at this build
	boxlo_old[0] = boxlo[0];
	boxlo_old[1] = boxlo[1];
	boxlo_old[2] = boxlo[2];
	boxhi_old[0] = boxhi[0];
	boxhi_old[1] = boxhi[1];
	boxhi_old[2] = boxhi[2];

	
	if (dist_check){
		// store coordinate at this build
		double **x = particle->x;
		for (int i = 0; i < nlocal; i++)
		for (int j = 0; j < 3; j++) {
			x_old[i][j] = x[i][j];
		}

		// GPU data transfer
		int n_master;
		if (slave_flag)
			n_master = particle->nlocal - slave_num;
		else
			n_master = particle->nlocal;
		gpuTransfer << < int(n_master + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> > (particle->devCoordXold, particle->devCoordX, n_master);
		gpuTransfer << < int(n_master + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> > (particle->devCoordYold, particle->devCoordY, n_master);
		gpuTransfer << < int(n_master + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> > (particle->devCoordZold, particle->devCoordZ, n_master);

	}
	build_before = 1;
	timer->stamp();
	create_linked_list();
	timer->stamp(TIME_NEIGHBOR_L);
    create_neigh_list();
	timer->stamp(TIME_NEIGHBOR_B);

}

/* ---------------------------------------------------------------------- 
   Allocate for lists
---------------------------------------------------------------------- */

void Neighbor::allocate(int nmax)
{
	// allocate NeighList class

	neighlist->grow(nmax);
	int nlocal = particle->nlocal;
	if (dist_check) {
		if (nlocal > nlocal_old) {
			nlocal_old = nlocal;
			memory->destroy(x_old);
			memory->create(x_old, nlocal + 10000, 3, "Neighbor: x_old");
		}
	}

	if (nmax >= nmax_old) {
		nmax_old = nmax;
		memory->destroy(linked_list);
		memory->create(linked_list, nmax_old, "Neighbor: linked_list");
	}
}

/* ---------------------------------------------------------------------- 
   Create linked list
   (Check the exact algorithm in the PDPS manual)
---------------------------------------------------------------------- */

void Neighbor::create_linked_list()
{
	int c;
	double *devCoordX = particle->devCoordX;
	double *devCoordY = particle->devCoordY;
	double *devCoordZ = particle->devCoordZ;
	double **x = particle->x;
	int *mask = particle->mask;
	int nall = particle->nlocal + particle->nghost;

	// Reset the headers to -1
	if (slave_flag){
		if (slave_linked){
			for (c = 0; c < subncxyz; c++) {
				head[c] = head_slave[c];
			}
		}
		else{
			for (c = 0; c < subncxyz; c++) {
				head_slave[c] = -1;
				head[c] = -1;
			}
			for (int i = nall - 1; i >= 0; i--) {
				if (mask[i] & slave_bit){
					c = coord2cell(x[i]);
					if (c < 0 || c >= subncxyz) {
						printf("c = %f", c);
						error->all(FLERR, "Invalid cell setup");
					}
					// Link to the previoius occupant (or -1 if you are the 1st)
					linked_list[i] = head_slave[c];
					// The last one goes to the header
					head[c] = head_slave[c] = i;
				}	//  mask[i] & slave_bit
			}	//  i = nall -1
			slave_linked = 1;
		}	// slave_linked
	}
	else{
		for (c = 0; c < subncxyz; c++) {
			head[c] = -1;
		}
	}		//  slave_flag

	// Scan particles to construct headers, head, & linked lists, linked_list
	// Scan ghost particles first, so they will appear at the end of the neighbor list
	// Scan in an reverse order, so they will apear in a correct order in the neighbor list 
	if (slave_flag){
		for (int i = nall - 1; i >= 0; i--) {
			if (!(mask[i] & slave_bit)){
				if (i == 0)
					i = i;
				c = coord2cell(x[i]);
				if (c < 0 || c >= subncxyz) {
					printf("c = %f", c);
					error->all(FLERR, "Invalid cell setup");
				}
				// Link to the previoius occupant (or -1 if you are the 1st)
				linked_list[i] = head[c];
				// The last one goes to the header
				head[c] = i;
			}		//  ! (mask[i] & slave_bit[i])

		}	//  int i

	}
	else{
		for (int i = nall - 1; i >= 0; i--) {
			c = coord2cell(x[i]);
			if (c < 0 || c >= subncxyz) {
				printf("c = %f", c);
				error->all(FLERR, "Invalid cell setup");
			}
			// Link to the previoius occupant (or -1 if you are the 1st)
			linked_list[i] = head[c];
			// The last one goes to the header
			head[c] = i;
		}
	}
	for (int i = nall - 1; i >= 0; i--){
		if (linked_list[i] < -1){
			c = coord2cell(x[i]);
			i = i;
		}
			
	}
	
	cudaError_t error_t;
	error_t = cudaMemcpy(devLinked_list, linked_list, (particle->nlocal + particle->nghost) * sizeof(int), cudaMemcpyHostToDevice);
	error_t = cudaMemcpy(devHead, head, subncxyz * sizeof(int), cudaMemcpyHostToDevice);

	/*gpulinklist << < int(nall + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >> > (devLinked_list, devHead, devSubclo, devSubnc,
		devBoxhi, devBoxlo, devCle, devNc, nall, subncxyz, subncxy, devCoordX, devCoordY, devCoordZ);*/

}

/* ----------------------------------------------------------------------
   map particles' coords into local cell index
   for orthogonal, only ghost atoms will have coord >= bboxhi or coord < bboxlo
     take special care to insure ghosts are in correct bins even w/ roundoff
     hi ghost atoms = nbin,nbin+1,etc
     owned atoms = 0 to nbin-1
     lo ghost atoms = -1,-2,etc
     this is necessary so that both procs on either side of PBC
       treat a pair of atoms straddling the PBC in a consistent way
   for triclinic, doesn't matter since stencil & neigh list built differently
------------------------------------------------------------------------- */

int Neighbor::coord2cell(double *x)
{
	int ic[3];

	for (int i = 0; i < 3; i++) {
		if (x[i] >= boxhi[i])
			ic[i] = static_cast<int> ((x[i] - boxhi[i])/cle[i]) + nc[i];
		else if (x[i] >= boxlo[i]) {
			ic[i] = static_cast<int> ((x[i] - boxlo[i])/cle[i]);
			ic[i] = MIN(ic[i],nc[i]-1);
		} else
			ic[i] = static_cast<int> ((x[i] - boxlo[i])/cle[i]) - 1;
	}

	if (domain->dim == 2)
		ic[2] = 0;

	for (int i = 0; i < 3; i++) {
		ic[i] -= subclo[i];
	}

	int c = (ic[2] + 2)*subncxy + (ic[1] + 2) * subnc[0] + ic[0] + 2;
	return c;
}

/* ----------------------------------------------------------------------
   Same as above, but return the cell location index
------------------------------------------------------------------------- */

int Neighbor::coord2cell(double *x, int &icx, int &icy, int &icz)
{
	int ic[3];

	for (int i = 0; i < 3; i++) {
		if (x[i] >= boxle[i])
			ic[i] = static_cast<int> ((x[i] - boxle[i])/cle[i]) + nc[i];
		else if (x[i] >= boxlo[i]) {
			ic[i] = static_cast<int> ((x[i] - boxlo[i])/cle[i]);
			ic[i] = MIN(ic[i],nc[i]-1);
		} else
			ic[i] = static_cast<int> ((x[i] - boxlo[i])/cle[i]) - 1;
	}

	if (domain->dim == 2)
		ic[2] = 0;

	for (int i = 0; i < 3; i++) {
		ic[i] -= subclo[i];
	}

	icx = ic[0];
	icy = ic[1];
	icz = ic[2];

	int c = ic[2]*subncxy + ic[1]*subnc[0] + ic[0];
	return c;
}

/* ---------------------------------------------------------------------- 
   Create neighbor list for all particles
   (Check the exact algorithm in the PDPS manual)
---------------------------------------------------------------------- */

void Neighbor::create_neigh_list()
{
	
	cudaError_t  error_t;
	

	int nall = particle->nlocal + particle->nghost;
	//(this->*pair_build)(neighlist);
	//error_t = cudaMemcpy(hostCoordX, particle->devCoordX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(hostCoordY, particle->devCoordY, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(hostCoordZ, particle->devCoordZ, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	//error_t = cudaMemcpy(hostNumneigh, devNumneigh, particle->nlocal * sizeof(int), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(hostPairtable, devPairtable, particle->nlocal * NEIGHMAX * sizeof(int), cudaMemcpyDeviceToHost);

	if (slave_flag){
		cudaEventRecord(start, 0);
		gpuneighbuild << < int(nall + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >(devLinked_list, devHead, devNumneigh, devPairtable, devCoffsets, noffsets,
			devSubclo, devSubnc, devBoxhi, devBoxlo, devCle, devNc, devRneighsq, particle->devMask,
			nall, subncxyz, subncxy, particle->devCoordX, particle->devCoordY, particle->devCoordZ, particle->devType, slave_bit);

	
		/*gpuneighbuild2 << < int(nall + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >> >(devNumneigh, devPairtable, devRneighsq, rneigh_max, particle->devMask,
			nall, particle->devCoordX, particle->devCoordY, particle->devCoordZ, particle->devType, slave_bit);*/
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time, start, stop);
	}
	else{
		gpuneighbuild << < int(nall + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >(devLinked_list, devHead, devNumneigh, devPairtable, devCoffsets, noffsets,
			devSubclo, devSubnc, devBoxhi, devBoxlo, devCle, devNc, devRneighsq, particle->devMask,
			nall, subncxyz, subncxy, particle->devCoordX, particle->devCoordY, particle->devCoordZ, particle->devType, 0);
	/*	gpuneighbuild2 << < int(nall + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >(devNumneigh, devPairtable, devRneighsq, rneigh_max, particle->devMask,
			nall, particle->devCoordX, particle->devCoordY, particle->devCoordZ, particle->devType, 0);*/

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time, start, stop);
	}

	//error_t = cudaMemcpy(hostNumneigh, devNumneigh, particle->nlocal * sizeof(int), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(hostPairtable, devPairtable, particle->nlocal * NEIGHMAX * sizeof(int), cudaMemcpyDeviceToHost);
	
	
}

/* ---------------------------------------------------------------------- 
   Check cell to see if it satisfies the boundary condtions
		: return 0 if the cell index is out of boundary
		  return 1 if the cell index is within the boundary
---------------------------------------------------------------------- */

int Neighbor::check_cell_boundary(int cell[3])
{
	for (int i = 0; i < 3; i++) {
		if (cell[i] < 0 || cell[i] >= subnc[i]) {
			return 0;
		}
	}
	return 1;
}

/* ---------------------------------------------------------------------- 
   Decide whether the neighbor list needs to be updated or not 
---------------------------------------------------------------------- */

int Neighbor::decide()
{
	nflag = 0;
	
	ago++;

	// during the simulation
	if (ago >= delay && ago % every == 0) {
		if (dist_check == 0) {
			nflag = 1;
		}
		else {
			nflag = check_distance();
		}
	}

	return nflag;
}

/* ----------------------------------------------------------------------
   if any particle moved trigger distance (half of neighbor skin) return 1
   shrink trigger distance if box size has changed
   conservative shrink procedure:
     compute distance each of 8 corners of box has moved since last reneighbor
     reduce skin distance by sum of 2 largest of the 8 values
     new trigger = 1/2 of reduced skin distance
   for orthogonal box, only need 2 lo/hi corners
   for triclinic, need all 8 corners since deformations can displace all 8
------------------------------------------------------------------------- */

int Neighbor::check_distance()
{
	double delx, dely, delz, rsq;
	double delta, deltasq, delta1, delta2;

	if (boxcheck) {
		delx = boxlo[0] - boxlo_old[0];
		dely = boxlo[1] - boxlo_old[1];
		delz = boxlo[2] - boxlo_old[2];
		delta1 = sqrt(delx*delx + dely*dely + delz*delz);
		delx = boxhi[0] - boxhi_old[0];
		dely = boxhi[1] - boxhi_old[1];
		delz = boxhi[2] - boxhi_old[2];
		delta2 = sqrt(delx*delx + dely*dely + delz*delz);
		delta = 0.5 * (rskin - (delta1 + delta2));
		deltasq = delta*delta;
	}
	else deltasq = triggersq;

	double **x = particle->x;
	int nlocal = particle->nlocal;

	int flag = 0;
	
	//	flag in GPU
	int *build_flag;
	cudaMalloc(&build_flag, sizeof(int));
	cudaMemcpy(build_flag, &flag, sizeof(int), cudaMemcpyHostToDevice);

	//for (int i = 0; i < nlocal; i++) {
	//	delx = x[i][0] - x_old[i][0];
	//	dely = x[i][1] - x_old[i][1];
	//	delz = x[i][2] - x_old[i][2];
	//	rsq = delx*delx + dely*dely + delz*delz;
	//	if (rsq > deltasq) {
	//		flag = 1;
	//	}
	//}
	
	//int flagall;
	//MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,mworld);
	//if (flagall && ago == MAX(every,delay)) {
	//	ndanger++;
	//	ndanger_total++;
	//	if (procid == 0) {
	//		fprintf(stdout, "Warning: ntimestep = %d ndanger = %d\n", update->ntimestep, ndanger);
	//		fflush(stdout);
	//	}
	//	//error->all(FLERR,"Dangerous build\n");
	//}
	
	int n_master;
	if (slave_flag)
		n_master = particle->nlocal - slave_num;
	else
		n_master = particle->nlocal;
	//cudaEvent_t start, stop;
	//float time;
	//cudaEventCreate(&start);
	//cudaEventCreate(&stop);
	//cudaEventRecord(start, 0);
	//cudaError_t error_t;
	//error_t = cudaMemcpy(neighbor->hostCoordX, particle->devCoordXold, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	//error_t = cudaMemcpy(neighbor->hostCoordZ, particle->devCoordX, particle->nlocal * sizeof(double), cudaMemcpyDeviceToHost);
	gpucheck_distance << < int (n_master + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> > (particle->devCoordX, particle->devCoordY, particle->devCoordZ,
		particle->devCoordXold, particle->devCoordYold, particle->devCoordZold,
		n_master, deltasq, build_flag);
	cudaMemcpy(&flag, build_flag, sizeof(int), cudaMemcpyDeviceToHost);
	//cudaEventRecord(stop, 0);
	//cudaEventSynchronize(stop);
	//cudaEventElapsedTime(&time, start, stop);
	return flag;
}

/* ---------------------------------------------------------------------- 
   Find distance 
---------------------------------------------------------------------- */

double Neighbor::find_distance(int id1, int id2)
{
	double distsq;
	double x1[3], x2[3];

	distsq = 0.0;
	for (int i = 0; i < 3; i++) {
		x1[i] = particle->x[id1][i];
		x2[i] = particle->x[id2][i];
		distsq += pow((x2[i] - x1[i]),2);
	}

	return distsq;
}

/* ---------------------------------------------------------------------- 
   Neighbor list debug (Only for nprocs = 1)
---------------------------------------------------------------------- */

void Neighbor::debug()
{
	FILE *file1;
	file1 = fopen("neighbor_list.txt","w");

	fprintf(file1, "TIMESTEP = %d\n", update->ntimestep);
	fprintf(file1, "local_id global_id x y\n");

	for (int i = 0; i < particle->nlocal + particle->nghost; i++) {
		fprintf(file1, "%d %d %g %g\n",i, particle->tag[i]-1,particle->x[i][0],particle->x[i][1]);
	}
	fprintf(file1,"\n");
	fclose(file1);
}

/* ----------------------------------------------------------------------
		set slaves particles which have no neighbor list 
---------------------------------------------------------------------- */

void Neighbor::setslave(int narg, char **arg)
{
	if (narg != 1) error->all(FLERR, "Illegal neighbor command");
	int slaveid = group->find_group(arg[0]);
	slave_bit = group->bitmask[slaveid];
	slave_flag = 1;
	slave_num = group->gparticles[slaveid];
}

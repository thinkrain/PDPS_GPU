/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "compute.h"
#include "domain.h"
#include "dump.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "parallel.h"
#include "particle.h"
#include "particle_type.h"
#include "update.h"

using namespace PDPS_NS;

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000
#define BIG 1.0e20

enum{SINGLE,MULTI};

Parallel::Parallel(PDPS *ps) : Pointers(ps)
{
	MPI_Comm_size(mworld,&nprocs);
	MPI_Comm_rank(mworld,&procid);

	grid2proc = NULL;
	procfactors = NULL;
	xsplit = NULL;
	ysplit = NULL;
	zsplit = NULL;
	bordergroup = 0;

	// let the sytle of parallelism the same as the neighbor list

	user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;
	for (int i = 0; i < 3; i++) {
		procgrid[i] = 1;
		procloc[i] = 0;
	}

	firstrecv = NULL;
	maxsendlist = NULL;

	pbc = NULL;
	pbc_flag = NULL;
	size_forward_recv = NULL;
	size_reverse_send = NULL;
	size_reverse_recv = NULL;
	slablo = NULL;
	slabhi = NULL; 

	multilo = multihi = NULL;
	rcghost_multi = NULL;

	buf_send = NULL;
	buf_recv = NULL;
	sendlist = NULL;
	maxsendlist = NULL;

	ghost_velocity = 1;    // default should be 0 

	comm_x_only = 1;
	comm_f_only = 1;

	maxsend = BUFMIN;
	memory->create(buf_send,maxsend+BUFEXTRA,"parallel: buf_send");
	maxrecv = BUFMIN;
	memory->create(buf_recv,maxrecv,"parallel: buf_recv");

	maxswap = 6;
	allocate_swap(maxswap);
	
	sendlist = (int **) memory->smalloc(maxswap*sizeof(int *),"Parallel: sendlist");
	memory->create(maxsendlist,maxswap,"Parallel: maxsendlist");
	for (int i = 0; i < maxswap; i++) {
		maxsendlist[i] = BUFMIN;
		sendlist[i] = NULL;
		memory->create(sendlist[i],BUFMIN,"comm:sendlist[i]");
	}
}

/* ---------------------------------------------------------------------- */

Parallel::~Parallel()
{
	memory->destroy(xsplit);
    memory->destroy(ysplit);
    memory->destroy(zsplit);
	
	memory->destroy(grid2proc);

	free_swap();

	if (sendlist) {
		for (int i = 0; i < maxswap; i++)
			memory->destroy(sendlist[i]);
	}

	memory->destroy(buf_send);
	memory->destroy(buf_recv);

	memory->sfree(sendlist);
	memory->destroy(maxsendlist);
}

/* ---------------------------------------------------------------------- */

void Parallel::init()
{
	style = neighbor->style; 
	map_style = particle->map_style;

	// comm_only = 1 if only x,f are exchanged in forward/reverse comm
    // comm_x_only = 0 if ghost_velocity since velocities are added

	comm_x_only = particle->ptype->comm_x_only;
	comm_f_only = particle->ptype->comm_f_only;
    if (ghost_velocity) comm_x_only = 0;
	
	// set per-particle sizes for forward/reverse/border comm
	// augment by velocity quantities if needed

	size_forward = particle->ptype->size_forward;
	size_reverse = particle->ptype->size_reverse;
	size_border = particle->ptype->size_border;

	if (ghost_velocity) {
		size_forward += particle->ptype->size_velocity;
		size_border += particle->ptype->size_velocity;
	}

	// maxforward = # of datums in largest forward communication
	// maxreverse = # of datums in largest reverse communication
	// query pair,fix,compute,dump for their requirements
	// pair style can force reverse comm even if newton off

	maxforward = MAX(size_forward,size_border);
	maxreverse = size_reverse;

	for (int i = 0; i < force->npairs; i++) {
		maxforward = MAX(maxforward,force->pair[i]->comm_forward);
		maxreverse = MAX(maxreverse,force->pair[i]->comm_reverse);
	}

	for (int i = 0; i < modify->nfixes; i++) {
		maxforward = MAX(maxforward,modify->fix[i]->comm_forward);
		maxreverse = MAX(maxreverse,modify->fix[i]->comm_reverse);
	}

	for (int i = 0; i < modify->ncomputes; i++) {
		maxforward = MAX(maxforward,modify->compute[i]->comm_forward);
		maxreverse = MAX(maxreverse,modify->compute[i]->comm_reverse);
	}

	for (int i = 0; i < output->ndumps; i++) {
		maxforward = MAX(maxforward,output->dump[i]->comm_forward);
		maxreverse = MAX(maxreverse,output->dump[i]->comm_reverse);
	}

	// memory for multi-style communication
	if (style == MULTI && multilo == NULL) {
		allocate_multi(maxswap);
		memory->create(rcghost_multi,particle->ntypes+1,3,"parallel: rcghost_multi");
	}
	if (style == SINGLE && multilo) {
		free_multi();
		memory->destroy(rcghost_multi);
	}
}

/* ----------------------------------------------------------------------
   Setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size
   single style sets slab boundaries (slablo,slabhi) based on max cutoff
   multi style sets type-dependent slab boundaries (multilo,multihi)
------------------------------------------------------------------------- */

void Parallel::setup()
{
	// cutghost[] = max distance at which ghost particles need to be acquired
	// for orthogonal:
	//   cutghost is in box coords = neigh->cutghost in all 3 dims
	// for triclinic:
	//   neigh->cutghost = distance between tilted planes in box coords
	//   cutghost is in lamda coords = distance between those planes
	// for multi:
	//   cutghostmulti = same as cutghost, only for each particle type

	int i;
	int ntypes = particle->ntypes;
	double *boxle, *sublo, *subhi;

	double rneigh;

	rneigh = neighbor->rneigh_max;
	if (rneigh == 0.0)
		error->all(FLERR,"the max cut-off has not been setup yet\n");

	boxle = domain->boxle;
	sublo = domain->sublo;
	subhi = domain->subhi;

	rcghost[0] = rcghost[1] = rcghost[2] = rneigh;
	
	if (style == MULTI) {
		for (i = 1; i <= ntypes; i++) {
			rcghost_multi[i][0] = rcghost_multi[i][1] =
				rcghost_multi[i][2] = neighbor->rneigh_type[i];
		}
	}

	// recvneed[idim][0/1] = # of procs away I recv particles from, within cutghost
	//   0 = from left, 1 = from right
	//   do not cross non-periodic boundaries, need[2] = 0 for 2d
	// sendneed[idim][0/1] = # of procs away I send particles to
	//   0 = to left, 1 = to right
	//   set equal to recvneed[idim][1/0] of neighbor proc
	// maxneed[idim] = max procs away any proc recvs particles in either direction
	// uniform = 1 = uniform sized sub-domains:
	//   maxneed is directly computable from sub-domain size
	//     limit to procgrid-1 for non-PBC
	//   recvneed = maxneed except for procs near non-PBC
	//   sendneed = recvneed of neighbor on each side
	// uniform = 0 = non-uniform sized sub-domains:
	//   compute recvneed via updown() which accounts for non-PBC
	//   sendneed = recvneed of neighbor on each side
	//   maxneed via Allreduce() of recvneed

	int *periodicity = domain->periodicity;
	int left, right;

	maxneed[0] = static_cast<int> (rcghost[0] * procgrid[0] / boxle[0]) + 1;
    maxneed[1] = static_cast<int> (rcghost[1] * procgrid[1] / boxle[1]) + 1;
    maxneed[2] = static_cast<int> (rcghost[2] * procgrid[2] / boxle[2]) + 1;

	if (domain->dim == 2) maxneed[2] = 0;

	// If rs is very large and cross the non periodic boundary 
	if (!periodicity[0]) maxneed[0] = MIN(maxneed[0],procgrid[0]-1);
    if (!periodicity[1]) maxneed[1] = MIN(maxneed[1],procgrid[1]-1);
    if (!periodicity[2]) maxneed[2] = MIN(maxneed[2],procgrid[2]-1);

	if (!periodicity[0]) {
      recvneed[0][0] = MIN(maxneed[0],procloc[0]);
      recvneed[0][1] = MIN(maxneed[0],procgrid[0]-procloc[0]-1);
      left = procloc[0] - 1;
      if (left < 0) left = procgrid[0] - 1;
      sendneed[0][0] = MIN(maxneed[0],procgrid[0]-left-1);
      right = procloc[0] + 1;
      if (right == procgrid[0]) right = 0;
      sendneed[0][1] = MIN(maxneed[0],right);
    } 
	else 
		recvneed[0][0] = recvneed[0][1] =
             sendneed[0][0] = sendneed[0][1] = maxneed[0];

    if (!periodicity[1]) {
      recvneed[1][0] = MIN(maxneed[1],procloc[1]);
      recvneed[1][1] = MIN(maxneed[1],procgrid[1]-procloc[1]-1);
      left = procloc[1] - 1;
      if (left < 0) left = procgrid[1] - 1;
      sendneed[1][0] = MIN(maxneed[1],procgrid[1]-left-1);
      right = procloc[1] + 1;
      if (right == procgrid[1]) right = 0;
      sendneed[1][1] = MIN(maxneed[1],right);
    } 
	else 
		recvneed[1][0] = recvneed[1][1] =
             sendneed[1][0] = sendneed[1][1] = maxneed[1];

    if (!periodicity[2]) {
      recvneed[2][0] = MIN(maxneed[2],procloc[2]);
      recvneed[2][1] = MIN(maxneed[2],procgrid[2]-procloc[2]-1);
      left = procloc[2] - 1;
      if (left < 0) left = procgrid[2] - 1;
      sendneed[2][0] = MIN(maxneed[2],procgrid[2]-left-1);
      right = procloc[2] + 1;
      if (right == procgrid[2]) right = 0;
      sendneed[2][1] = MIN(maxneed[2],right);
    } 
	else 
		recvneed[2][0] = recvneed[2][1] =
             sendneed[2][0] = sendneed[2][1] = maxneed[2];

	// allocate parallel memory
	
	nswap = 2 * (maxneed[0] + maxneed[1] + maxneed[2]);
	if (nswap > maxswap) grow_swap(nswap);

	// setup parameters for each exchange:
	// sendproc = proc to send to at each swap
	// recvproc = proc to recv from at each swap
	// for style SINGLE:
	//   slablo/slabhi = boundaries for slab of particles to send at each swap
	//   use -BIG/midpt/BIG to insure all particles included even if round-off occurs
	//   if round-off, particles recvd across PBC can be < or > than subbox boundary
	//   note that borders() only loops over subset of particles during each swap
	//   treat all as PBC here, non-PBC is handled in borders() via r/s need[][]
	// for style MULTI:
	//   multilo/multihi is same, with slablo/slabhi for each particle type
	// pbc_flag: 0 = nothing across a boundary, 1 = something across a boundary
	// pbc = -1/0/1 for PBC factor in each of 3/6 orthogonal/triclinic dirs
	// for triclinic, slablo/hi and pbc_border will be used in lamda (0-1) coords
	// 1st part of if statement is sending to the west/south/down
	// 2nd part of if statement is sending to the east/north/up

	int dim,ineed;
		
	int iswap = 0;
	for (dim = 0; dim < 3; dim++) {
		for (ineed = 0; ineed < 2*maxneed[dim]; ineed++) {
			pbc_flag[iswap] = 0;
			pbc[iswap][0] = pbc[iswap][1] = pbc[iswap][2] =
			  pbc[iswap][3] = pbc[iswap][4] = pbc[iswap][5] = 0;

			if (ineed % 2 == 0) {
				sendproc[iswap] = procneigh[dim][0];
				recvproc[iswap] = procneigh[dim][1];
				if (style == SINGLE) {
					if (ineed < 2) slablo[iswap] = -BIG;
					else slablo[iswap] = 0.5 * (sublo[dim] + subhi[dim]);
					slabhi[iswap] = sublo[dim] + rcghost[dim];
				}
				else {
					for (i = 1; i <= ntypes; i++) {
						if (ineed < 2) multilo[iswap][i] = -BIG;
						else multilo[iswap][i] = 0.5 * (sublo[dim] + subhi[dim]);
						multihi[iswap][i] = sublo[dim] + rcghost_multi[i][dim];
					 }
				}
				if (procloc[dim] == 0) {
					pbc_flag[iswap] = 1;
					pbc[iswap][dim] = 1;
				}
			} // if (ineed % 2 == 0)
			else { 
				sendproc[iswap] = procneigh[dim][1];
				recvproc[iswap] = procneigh[dim][0];
				if (style == SINGLE) {
					slablo[iswap] = subhi[dim] - rcghost[dim];
					if (ineed < 2) slabhi[iswap] = BIG;
					else slabhi[iswap] = 0.5 * (sublo[dim] + subhi[dim]);
				}
				else {
					for (i = 1; i <= ntypes; i++) {
						multilo[iswap][i] = subhi[dim] - rcghost_multi[i][dim];
						if (ineed < 2) multihi[iswap][i] = BIG;
						else multihi[iswap][i] = 0.5 * (sublo[dim] + subhi[dim]);
					 }
				} 
				if (procloc[dim] == (procgrid[dim] - 1)) {
					pbc_flag[iswap] = 1;
					pbc[iswap][dim] = -1;
				}
			} // if (ineed % 2 != 0)
			iswap++;
		} // for (ineed = 0; ineed < 2*maxneed[dim]; ineed++)
	} // for (dim = 0; dim < 3; dim++)
}

/* ---------------------------------------------------------------------- */

void Parallel::set_processors(int narg, char **arg)
{
	if (narg < 3) error->all(FLERR, "Illegal processors command");

	if (strcmp(arg[0], "*") == 0) user_procgrid[0] = 0;
	else user_procgrid[0] = atoi(arg[0]);
	if (strcmp(arg[1], "*") == 0) user_procgrid[1] = 0;
	else user_procgrid[1] = atoi(arg[1]);
	if (strcmp(arg[2], "*") == 0) {
		if (domain->dim == 2) user_procgrid[2] = 1;
		else user_procgrid[2] = 0;
	}
	else user_procgrid[2] = atoi(arg[2]);

	if (user_procgrid[0] < 0 || user_procgrid[1] < 0 || user_procgrid[2] < 0){
		error->all(FLERR,"Illegal processors command");
	}

	if (domain->dim == 2 && user_procgrid[2] > 1) {
		error->all(FLERR, "Only one processor is allowed for 2D problem");
	}

	int p = user_procgrid[0]*user_procgrid[1]*user_procgrid[2];
	if (p > 0 && p != nprocs) {
		error->all(FLERR,"User specified processors != physical requested processors");
	}
}

/* ---------------------------------------------------------------------- */

void Parallel::set_proc_grid(int outflag)
{
	int npossible;
	int dim = domain->dim;

	if (dim == 3) {
		npossible = find_factors3(nprocs,NULL);
		memory->create(procfactors,npossible,3,"Parallel: procfactors");
		npossible = find_factors3(nprocs,procfactors);
	}
	if (dim == 2) {
		npossible = find_factors2(nprocs,NULL);
		memory->create(procfactors,npossible,3,"Parallel: procfactors");
		npossible = find_factors2(nprocs,procfactors);
	}

	npossible = filter_user(npossible, procfactors, 3, user_procgrid);

	if (npossible == 0) {
		error->all(FLERR,"Could not create 3d grid of processors");
	}

	best_factors(npossible, procfactors, procgrid, 1, 1, 1);

	

	// clearn-up

	memory->destroy(procfactors);

	if (procgrid[0]*procgrid[1]*procgrid[2] != nprocs) {
		char str[128];
		sprintf(str,"Invalid grid of %d processors",nprocs);
		error->all(FLERR,str);
	}

	// call MPI_Cart to gather information of processor grid	
	if (grid2proc) memory->destroy(grid2proc);
	memory->create(grid2proc,procgrid[0],procgrid[1],procgrid[2],
                 "comm:grid2proc");

	cart_map(0,procgrid,procloc,procneigh,grid2proc);

	// print 3d grid of processors information to screen and logfile

	if (procid == 0) {
		if (screen) {
			fprintf(screen,"  %d by %d by %d MPI processor grid\n",
              procgrid[0],procgrid[1],procgrid[2]);
		}
		if (logfile) {
			fprintf(logfile,"  %d by %d by %d MPI processor grid\n",
              procgrid[0],procgrid[1],procgrid[2]);
		}
	}

	// set xplist, ysplit and zsplit for uniform spacings

	memory->destroy(xsplit);
	memory->destroy(ysplit);
	memory->destroy(zsplit);

	memory->create(xsplit,procgrid[0]+1,"Parallel: xsplit");
	memory->create(ysplit,procgrid[1]+1,"Parallel: ysplit");
	memory->create(zsplit,procgrid[2]+1,"Parallel: zsplit");

	for (int i = 0; i < procgrid[0]; i++) xsplit[i] = i * 1.0/procgrid[0];
	for (int i = 0; i < procgrid[1]; i++) ysplit[i] = i * 1.0/procgrid[1];
	for (int i = 0; i < procgrid[2]; i++) zsplit[i] = i * 1.0/procgrid[2];

	xsplit[procgrid[0]] = ysplit[procgrid[1]] = zsplit[procgrid[2]] = 1.0;
}

/* ----------------------------------------------------------------------
   generate all possible 3-integer factorizations of N
   store them in factors if non-NULL
   return # of factorizations
------------------------------------------------------------------------- */

int Parallel::find_factors3(int n, int **factors)
{
	int i, j, nyz;

	int m = 0;

	for (i = 1; i <= n; i++) {
		if (n % i != 0) continue;
		nyz = n/i;
		for (j = 1; j <= nyz; j++) {
			if (nyz % j != 0) continue;
			if (factors != NULL) {
				factors[m][0] = i;
				factors[m][1] = j;
				factors[m][2] = nyz/j;
			}
			m++;
		}
	}

	return m;
}

/* ----------------------------------------------------------------------
   generate all possible 2-integer factorizations of N
   store them in factors if non-NULL
   return # of factorizations
------------------------------------------------------------------------- */

int Parallel::find_factors2(int n, int **factors)
{
	int i;

	int m = 0;

	for (i = 1; i <= n; i++) {
		if (n % i != 0) continue;

		if (factors != NULL) {
			factors[m][0] = i;
			factors[m][1] = n/i;
			factors[m][2] = 1;
		}
		m++;
	}

	return m;
}

/* ----------------------------------------------------------------------
   choose best factors from list of Npossible factors
   best = minimal surface area of sub-domain
   return best = 3 factors
   return index of best factors in factors
------------------------------------------------------------------------- */

int Parallel::best_factors(int npossible, int **factors, int *best,
                          const int sx, const int sy, const int sz)
{
	// determine cross-sectional areas for orthogonal and triclinic boxes
	// for triclinic, area = cross product of 2 edge vectors stored in h matrix
	// area[3] = surface area 3 box faces divided by sx,sy,sz
	// area[0] = xy, area[1] = xz, area[2] = yz

	double area[3];
	
	area[0] = domain->boxle[0] * domain->boxle[1] / (sx*sy);
	area[1] = domain->boxle[0] * domain->boxle[2] / (sx*sz);
	area[2] = domain->boxle[1] * domain->boxle[2] / (sy*sz);

	int index;
	double surf;
	double bestsurf = 2.0 * (area[0]+area[1]+area[2]);

	for (int m = 0; m < npossible; m++) {
		surf = area[0]/factors[m][0]/factors[m][1] +
		  area[1]/factors[m][0]/factors[m][2] +
		  area[2]/factors[m][1]/factors[m][2];
		if (surf < bestsurf) {
		  bestsurf = surf;
		  best[0] = factors[m][0];
		  best[1] = factors[m][1];
		  best[2] = factors[m][2];
		  index = m;
		}
	}

	return index;
}

/* ----------------------------------------------------------------------
   remove any factors that do not match non-zero user_factors Px,Py,Pz
------------------------------------------------------------------------- */

int Parallel::filter_user(int n, int **factors, int m, int *user_factors)
{
	int i = 0;
	while (i < n) {
		int flag = 0;
		if (user_factors[0] && factors[i][0] != user_factors[0]) flag = 1;
		if (user_factors[1] && factors[i][1] != user_factors[1]) flag = 1;
		if (user_factors[2] && factors[i][2] != user_factors[2]) flag = 1;
		if (flag) {
			for (int j = 0; j < m; j++) factors[i][j] = factors[n-1][j];
			n--;
		} else i++;
	}
	return n;
}

/* ----------------------------------------------------------------------
   map processors to 3d grid via MPI_Cart routines
   MPI may do layout in machine-optimized fashion
------------------------------------------------------------------------- */

void Parallel::cart_map(int reorder, int *procgrid,
                       int *procloc, int procneigh[3][2], int ***grid2procroc)
{
	int periods[3];
	periods[0] = periods[1] = periods[2] = 1;
	MPI_Comm comm_cart;

	MPI_Cart_create(mworld,3,procgrid,periods,reorder,&comm_cart);
	MPI_Cart_get(comm_cart,3,procgrid,periods,procloc);
	MPI_Cart_shift(comm_cart,0,1,&procneigh[0][0],&procneigh[0][1]);
	MPI_Cart_shift(comm_cart,1,1,&procneigh[1][0],&procneigh[1][1]);
	MPI_Cart_shift(comm_cart,2,1,&procneigh[2][0],&procneigh[2][1]);

	int coords[3];
	int i,j,k;
	for (i = 0; i < procgrid[0]; i++)
	for (j = 0; j < procgrid[1]; j++)
	  for (k = 0; k < procgrid[2]; k++) {
		coords[0] = i; coords[1] = j; coords[2] = k;
		MPI_Cart_rank(comm_cart,coords,&grid2procroc[i][j][k]);
	  }

	MPI_Comm_free(&comm_cart);
}

/* ----------------------------------------------------------------------
   exchange: move particles to correct processors
   particles exchanged with all 6 stencil neighbors
   send out particles that have left my box, receive ones entering my box
   particles will be lost if not inside some proc's box
     can happen if particle moves outside of non-periodic bounary
     or if particle moves more than one proc away
   this routine called before every reneighboring
   for triclinic, particles must be in lamda coords (0-1) before exchange is called
------------------------------------------------------------------------- */

void Parallel::exchange()
{
	int i, m, nsend, nrecv, nrecv1, nrecv2, nlocal;
	double lo, hi, value;
	double **x;
	double *sublo, *subhi, *buf;
	ParticleType *ptype = particle->ptype;

	MPI_Request request;
	MPI_Status status;

	// clear global->local map for owned and ghost particles
    // b/c particles migrate to new procs in exchange() and
    // new ghosts are created in borders()
    // map_set() is done at end of borders()
	if (map_style) particle->map_clear();

	sublo = domain->sublo;
	subhi = domain->subhi;

	// loop over dimensions

	for (int dim = 0; dim < 3; dim++) {

		// fill buffer with particles leaving my box, using < and >=
		// when particle is deleted, fill it in with last particle
	
		x = particle->x;
		lo = sublo[dim];
		hi = subhi[dim];
		nlocal = particle->nlocal;
		i = nsend = 0;
		// loop for local particles
		while (i < nlocal) {
			if (x[i][dim] < lo || x[i][dim] >= hi) {
				if (nsend > maxsend) grow_send(nsend,1);
				nsend += ptype->pack_exchange(i,&buf_send[nsend]);
				ptype->copyI2J(nlocal-1,i,1);
				nlocal--;
			} 
			else i++;
		}
		particle->nlocal = nlocal;
		// send/recv particles in both directions
		// if 1 proc in dimension, no send/recv, set recv buf to send buf
		// if 2 procs in dimension, single send/recv
		// if more than 2 procs in dimension, send/recv to both neighbors
		
		if (procgrid[dim] == 1) {
			nrecv = nsend;
			buf = buf_send;
		}
	
		else {
			MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,
				         &nrecv1,1,MPI_INT,procneigh[dim][1],0,mworld,&status);
			nrecv = nrecv1;
			if (procgrid[dim] > 2) {
				MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,
					         &nrecv2,1,MPI_INT,procneigh[dim][0],0,mworld,&status);
				nrecv += nrecv2;
			}
			if (nrecv > maxrecv) grow_recv(nrecv);
			MPI_Irecv(buf_recv,nrecv1,MPI_DOUBLE,procneigh[dim][1],0,
                mworld,&request);
			MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][0],0,mworld);
			MPI_Wait(&request,&status);
				
			if (procgrid[dim] > 2) {
				MPI_Irecv(&buf_recv[nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,
					  mworld,&request);
				MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][1],0,mworld);
				MPI_Wait(&request,&status);
			}

			buf = buf_recv;
		} // if (procgrid[dim] != 1)
		// check incoming particles to see if they are in my box
		// if so, add to my list
		m = 0;
		while (m < nrecv) {
			value = buf[m+dim+1];
			if (value >= lo && value < hi) 
				m += ptype->unpack_exchange(&buf[m]);
			else 
				m += static_cast<int> (buf[m]);
		}
	} // for (int dim = 0; dim < 3; dim++)
}

/* ----------------------------------------------------------------------
   borders: list nearby particles to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a communicate (so don't need to explicitly
     call communicate routine on reneighboring timestep)
   this routine is called before every reneighboring
------------------------------------------------------------------------- */

void Parallel::borders()
{
	int i, n, itype, iswap, dim, ineed, twoneed, smax, rmax;
	int nsend, nrecv, sendflag, nfirst, nlast, ngroup;
	double lo, hi;
	int *type;
	double **x;
	double *buf, *mlo, *mhi;
	MPI_Request request;
	MPI_Status status;
	ParticleType *ptype = particle->ptype;

	// clear old ghosts and any ghost bonus data internal to particleVec

	particle->nghost = 0;

	iswap = 0;
	smax = rmax = 0;

	for (dim = 0; dim < 3; dim++) {
		nlast = 0;
		twoneed = 2*maxneed[dim];
		for (ineed = 0; ineed < twoneed; ineed++) {
			
			// find particles within slab boundaries lo/hi using <= and >=
			// check particles between nfirst and nlast
			//   for first swaps in a dim, check owned and ghost
			//   for later swaps in a dim, only check newly arrived ghosts
			// store sent particle indices in list for use in future timesteps

			x = particle->x;
			if (style == SINGLE) {
				lo = slablo[iswap];
				hi = slabhi[iswap];
			}
			else {
				type = particle->type;
				mlo = multilo[iswap];
				mhi = multihi[iswap];
			}

			if (ineed % 2 == 0) {
				nfirst = nlast;
				nlast = particle->nlocal + particle->nghost;
			}

			nsend = 0;

			// sendflag = 0 if I do not send on this swap
			// sendneed test indicates receiver no longer requires data
			// e.g. due to non-PBC or non-uniform sub-domains

			if (ineed/2 >= sendneed[dim][ineed % 2]) sendflag = 0;
			else sendflag = 1;

			// find send particles according to SINGLE vs MULTI
			// all particles eligible versus particles in bordergroup
			// only need to limit loop to bordergroup for first sends (ineed < 2)
			// on these sends, break loop in two: owned (in group) and ghost

			if (sendflag) {
				if (!bordergroup || ineed >= 2) {
					if (style == SINGLE) {
						for (i = nfirst; i < nlast; i++) {
							if (x[i][dim] >= lo && x[i][dim] <= hi) {
								if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
								sendlist[iswap][nsend++] = i;
							}
						}
					}
					else {
						for (i = nfirst; i < nlast; i++) {
							itype = type[i];
							if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
								if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
								sendlist[iswap][nsend++] = i;
							}
						}
					}
				} // if (!bordergroup || ineed >= 2)
				else {
					error->all(FLERR,"This part has not been investigated yet");
					if (style == SINGLE) {
						ngroup = particle->nfirst; // ?
						for (i = 0; i < ngroup; i++) {
							if (x[i][dim] >= lo && x[i][dim] <= hi) {
								if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
								sendlist[iswap][nsend++] = i;
							}
						}
						for (i = particle->nlocal; i < nlast; i++) {
							if (x[i][dim] >= lo && x[i][dim] <= hi) {
								if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
								sendlist[iswap][nsend++] = i;
							}
						}
					}
					else {
						ngroup = particle->nfirst;
						for (i = 0; i < ngroup; i++) {
							itype = type[i];
							if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
								if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
								sendlist[iswap][nsend++] = i;
							}
						}
						for (i = particle->nlocal; i < nlast; i++) {
							itype = type[i];
							if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
								if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
								sendlist[iswap][nsend++] = i;
							}
						}
					}
				} // else (!bordergroup || ineed >= 2)
			} // if (sendflag)

			// pack up list of border particles

			if (nsend*size_border > maxsend) {
				grow_send(nsend*size_border,0);
			}
			if (ghost_velocity) {
				 n = ptype->pack_border_vel(nsend,sendlist[iswap],buf_send,pbc_flag[iswap],pbc[iswap]);
			}
			else {
				 n = ptype->pack_border(nsend,sendlist[iswap],buf_send,
                              pbc_flag[iswap],pbc[iswap]);
			}

			// swap particles with other proc
		    // no MPI calls except SendRecv if nsend/nrecv = 0
		    // put incoming ghosts at end of my particle arrays
		    // if swapping with self, simply copy, no messages
			
			
			if (sendproc[iswap] != procid) {
				MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
							 &nrecv,1,MPI_INT,recvproc[iswap],0,mworld,&status);
				if (nrecv*size_border > maxrecv) grow_recv(nrecv*size_border);
				
				if (nrecv) MPI_Irecv(buf_recv,nrecv*size_border,MPI_DOUBLE,
									 recvproc[iswap],0,mworld,&request);
				if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,mworld);
				if (nrecv) MPI_Wait(&request,&status);
				buf = buf_recv;
			} 
			else {
				nrecv = nsend;
				buf = buf_send;
		    }
			
			// unpack buffer

			if (ghost_velocity)
				ptype->unpack_border_vel(nrecv,particle->nlocal+particle->nghost,buf);
			else
				ptype->unpack_border(nrecv,particle->nlocal+particle->nghost,buf);

			// set all pointers & counters

			smax = MAX(smax,nsend);
			rmax = MAX(rmax,nrecv);
			sendnum[iswap] = nsend;
			recvnum[iswap] = nrecv;
			size_forward_recv[iswap] = nrecv*size_forward;
			size_reverse_send[iswap] = nrecv*size_reverse;
			size_reverse_recv[iswap] = nsend*size_reverse;
			firstrecv[iswap] = particle->nlocal + particle->nghost;
			particle->nghost += nrecv;
			iswap++;

		} // for (ineed = 0; ineed < twoneed; ineed++)
	} // for (dim = 0; dim < 3; dim++)

	// insure send/recv buffers are long enough for all forward & reverse comm

	int max = MAX(maxforward*smax,maxreverse*rmax);
	if (max > maxsend) grow_send(max,0);
	max = MAX(maxforward*rmax,maxreverse*smax);
	if (max > maxrecv) grow_recv(max);

	// reset global->local map
	if (map_style) particle->map_set();
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void Parallel::grow_send(int n, int flag)
{
	maxsend = static_cast<int> (BUFFACTOR * n);
	if (flag)
		memory->grow(buf_send,(maxsend+BUFEXTRA),"Parallel: buf_send");
	else {
		memory->destroy(buf_send);
		memory->create(buf_send,maxsend+BUFEXTRA,"Parallel, buf_send");
	}
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void Parallel::grow_recv(int n)
{
	maxrecv = static_cast<int> (BUFFACTOR * n);
	memory->destroy(buf_recv);
	memory->create(buf_recv,maxrecv,"parallel: buf_recv");
}

/* ----------------------------------------------------------------------
   allocation of swap info
------------------------------------------------------------------------- */

void Parallel::allocate_swap(int n)
{
	memory->create(sendnum,n,"Parallel: sendnum");
	memory->create(recvnum,n,"Parallel: recvnum");
	memory->create(sendproc,n,"Parallel: sendproc");
	memory->create(recvproc,n,"Parallel: recvproc");
	memory->create(size_forward_recv,n,"Parallel: size");
	memory->create(size_reverse_send,n,"Parallel: size");
	memory->create(size_reverse_recv,n,"Parallel: size");
	memory->create(slablo,n,"Parallel: slablo");
	memory->create(slabhi,n,"Parallel: slabhi");
	memory->create(firstrecv,n,"Parallel: firstrecv");
	memory->create(pbc_flag,n,"Parallel: pbc_flag");
	memory->create(pbc,n,6,"Parallel: pbc");
}

/* ----------------------------------------------------------------------
   allocation array
------------------------------------------------------------------------- */

void Parallel::allocate_multi(int n)
{
	multilo = memory->create(multilo, n, particle->ntypes+1, "parallel: multilo");
	multihi = memory->create(multihi, n, particle->ntypes+1, "parallel: multihi");
}

/* ----------------------------------------------------------------------
   free memory for swaps
------------------------------------------------------------------------- */

void Parallel::free_swap()
{
	memory->destroy(sendnum);
	memory->destroy(recvnum);
	memory->destroy(sendproc);
	memory->destroy(recvproc);
	memory->destroy(size_forward_recv);
	memory->destroy(size_reverse_send);
	memory->destroy(size_reverse_recv);
	memory->destroy(slablo);
	memory->destroy(slabhi);
	memory->destroy(firstrecv);
	memory->destroy(pbc_flag);
	memory->destroy(pbc);
}

/* ----------------------------------------------------------------------
   free memory for swaps
------------------------------------------------------------------------- */

void Parallel::free_multi()
{
	memory->destroy(multilo);
	memory->destroy(multihi);
}

/* ----------------------------------------------------------------------
   realloc the buffers needed for swaps
------------------------------------------------------------------------- */

void Parallel::grow_swap(int n)
{
	free_swap();
	allocate_swap(n);
	if (style == MULTI) {
		free_multi();
		allocate_multi(n);
	}

	sendlist = (int **)
		memory->srealloc(sendlist,n*sizeof(int *),"paralle: sendlist");
	memory->grow(maxsendlist,n,"parallel: maxsendlist");
	for (int i = maxswap; i < n; i++) {
		maxsendlist[i] = BUFMIN;
		sendlist[i] = NULL;
		memory->create(sendlist[i],BUFMIN,"parallel: sendlist[i]");
	}
	maxswap = n;
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void Parallel::grow_list(int iswap, int n)
{
	maxsendlist[iswap] = static_cast<int> (BUFFACTOR * n);
	memory->grow(sendlist[iswap],maxsendlist[iswap],"Parallel: sendlist[iswap]");
}

/* ----------------------------------------------------------------------
   forward communication of particle coords every timestep
   other per-particle attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void Parallel::forward_comm(int dummy)
{
	int n;
	MPI_Request request;
	MPI_Status status;
	double **x = particle->x;
	double *buf;
	ParticleType *ptype = particle->ptype;

	// exchange data with another proc
	// if other proc is self, just copy
	// if comm_x_only set, exchange or copy directly to x, don't unpack

	for (int iswap = 0; iswap < nswap; iswap++) {
		if (sendproc[iswap] != procid) {
			if (comm_x_only) {
				if (size_forward_recv[iswap]) buf = x[firstrecv[iswap]];
				else buf = NULL;
				if (size_forward_recv[iswap])
					MPI_Irecv(buf,size_forward_recv[iswap],MPI_DOUBLE,
						recvproc[iswap],0,mworld,&request);
				n = ptype->pack_comm(sendnum[iswap],sendlist[iswap],
						buf_send,pbc_flag[iswap],pbc[iswap]);
				if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,mworld);
				if (size_forward_recv[iswap]) MPI_Wait(&request,&status);
			} 
			else if (ghost_velocity) {
				if (size_forward_recv[iswap])
					MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,recvproc[iswap],0,mworld,&request);
				n = ptype->pack_comm_vel(sendnum[iswap],sendlist[iswap],buf_send,pbc_flag[iswap],pbc[iswap]);
				if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,mworld);
				if (size_forward_recv[iswap]) MPI_Wait(&request,&status);
				ptype->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_recv);
			}
			else {
				if (size_forward_recv[iswap])
					MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
						recvproc[iswap],0,mworld,&request);
				if (n) MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,mworld);
				if (size_forward_recv[iswap]) MPI_Wait(&request,&status);
				ptype->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_recv);
			}
		} // if (sendproc[iswap] != procid)
		else {
			if (comm_x_only) {
				if (sendnum[iswap])
					n = ptype->pack_comm(sendnum[iswap],sendlist[iswap],
                              x[firstrecv[iswap]],pbc_flag[iswap],
                              pbc[iswap]);
			}
			else if (ghost_velocity) {
				n = ptype->pack_comm_vel(sendnum[iswap],sendlist[iswap],
                                buf_send,pbc_flag[iswap],pbc[iswap]);
				ptype->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_send);
			}
			else {
				n = ptype->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);
				ptype->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_send);
			}
		} // if (sendproc[iswap] = procid)
	} // for (int iswap = 0; iswap < nswap; iswap++)
}

/* ----------------------------------------------------------------------
   Reverse communication of forces on particles every timestep
   other per-particle attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void Parallel::reverse_comm()
{
	int n;
	MPI_Request request;
	MPI_Status status;
	double **f = particle->f;
	double *buf;
	ParticleType *ptype = particle->ptype;

	// exchange data with another proc
	// if other proc is self, just copy
	// if comm_f_only set, exchange or copy directly from f, don't pack

	for (int iswap = nswap-1; iswap >= 0; iswap--) {
		if (sendproc[iswap] != procid) {
			if (comm_f_only) {
				if (size_reverse_recv[iswap])
					MPI_Irecv(buf_recv,size_reverse_recv[iswap],MPI_DOUBLE,
						sendproc[iswap],0,mworld,&request);
				if (size_reverse_send[iswap]) buf = f[firstrecv[iswap]];
				else buf = NULL;
				if (size_reverse_send[iswap])
					MPI_Send(buf,size_reverse_send[iswap],MPI_DOUBLE,recvproc[iswap],0,mworld);
				if (size_reverse_recv[iswap]) MPI_Wait(&request,&status);
			}
			else {
				if (size_reverse_recv[iswap])
					MPI_Irecv(buf_recv,size_reverse_recv[iswap],MPI_DOUBLE,
						sendproc[iswap],0,mworld,&request);
				n = ptype->pack_reverse(recvnum[iswap],firstrecv[iswap],buf_send);
				if (n) MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,mworld);
				if (size_reverse_recv[iswap]) MPI_Wait(&request,&status);
			}
			ptype->unpack_reverse(sendnum[iswap],sendlist[iswap],buf_recv);
		} // if (sendproc[iswap] != procid)
		else {
			if (comm_f_only) {
				if (sendnum[iswap]) 
					ptype->unpack_reverse(sendnum[iswap],sendlist[iswap],
						f[firstrecv[iswap]]);
			}
			else {
				n = ptype->pack_reverse(recvnum[iswap],firstrecv[iswap],buf_send);
				ptype->unpack_reverse(sendnum[iswap],sendlist[iswap],buf_send);
			} 
		} // if (sendproc[iswap] = procid)
	} // for (int iswap = nswap-1; iswap >= 0; iswap--)
}

/* ----------------------------------------------------------------------
Reverse communication of rho on particles every timestep
other per-particle attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void Parallel::reverse_comm_pair(Pair *pair)
{
	int n;
	MPI_Request request;
	MPI_Status status;
	double *buf;

	// exchange data with another proc
	// if other proc is self, just copy
	// if comm_f_only set, exchange or copy directly from f, don't pack
	for (int iswap = nswap - 1; iswap >= 0; iswap--) {
		if (sendproc[iswap] != procid) {

			if (sendnum[iswap])
				MPI_Irecv(buf_recv, pair->comm_reverse * sendnum[iswap], MPI_DOUBLE, sendproc[iswap], 0, mworld, &request);
			n = pair->pack_reverse_comm(recvnum[iswap], firstrecv[iswap], buf_send);
			if (n)
			{
				MPI_Send(buf_send, n, MPI_DOUBLE, recvproc[iswap], 0, mworld);
			}
			if (sendnum[iswap])
				MPI_Wait(&request, &status);
			pair->unpack_reverse_comm(sendnum[iswap], sendlist[iswap], buf_recv);
		} // if (sendproc[iswap] != procid)
		else{
			n = pair->pack_reverse_comm(recvnum[iswap], firstrecv[iswap], buf_send);
			pair->unpack_reverse_comm(sendnum[iswap], sendlist[iswap], buf_send);
		}

	} // for (int iswap = nswap-1; iswap >= 0; iswap--)
}

/* ----------------------------------------------------------------------
Forward communication of rho on particles every timestep
other per-particle attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void Parallel::forward_comm_pair(Pair *pair)
{
	int n;
	MPI_Request request;
	MPI_Status status;
	double *buf;

	// exchange data with another proc
	// if other proc is self, just copy
	// if comm_f_only set, exchange or copy directly from f, don't pack
	for (int iswap = 0; iswap < nswap; iswap++) {

		if (sendproc[iswap] != procid) {

			if (recvnum[iswap])
				MPI_Irecv(buf_recv, pair->comm_forward * recvnum[iswap], MPI_DOUBLE, recvproc[iswap], 0, mworld, &request);
			n = pair->pack_forward_comm(sendnum[iswap], sendlist[iswap], buf_send);
			if (n)
			{
				MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, mworld);
			}
			if (recvnum[iswap])
				MPI_Wait(&request, &status);
			pair->unpack_forward_comm(recvnum[iswap], firstrecv[iswap], buf_recv);
		} // if (sendproc[iswap] != procid)
		else{
			n = pair->pack_forward_comm(sendnum[iswap], sendlist[iswap], buf_send);
			pair->unpack_forward_comm(recvnum[iswap], firstrecv[iswap], buf_send);
		}

	} // for (int iswap = iswap; iswap >= 0; iswap++)
}

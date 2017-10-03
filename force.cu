/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

// system library
#include "math.h"
#include "stdlib.h"
#include "string.h"

// pdps library
#include "error.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "pair.h"
#include "style_pair.h"
#include "particle.h"
#include "update.h"
#include "parallel.h"
#include "timer.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
using namespace PDPS_NS;

#define DELTA 1
#define BIG 1.0e20

__global__ void gpuforce_clear(const int nall, double *devForceX, double *devForceY, double *devForceZ){
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < nall){
		devForceX[i] = 0.0;
		devForceY[i] = 0.0;
		devForceZ[i] = 0.0;
	}

}

/* ---------------------------------------------------------------------- */

Force::Force(PDPS *ps) : Pointers(ps)
{
	maxpairs = 0;
	npairs = 0;

	// initialization
	boltz = 0.0;
	hplanck = 0.0;
	mvv2e = 0.0;
	ftm2v = 0.0;
	mv2d = 0.0;
	nktv2p = 0.0;
	qqr2e = 0.0;
	qe2f = 0.0;
	vxmu2f = 0.0;
	xxt2kmu = 0.0;
	dielectric = 0.0;
	qqrd2e = 0.0;
	e_mass = 0.0;
	hhmrr2e = 0.0;
	mvh2r = 0.0;

	pair = NULL;
	type2pair = NULL;
}

/* ---------------------------------------------------------------------- */

Force::~Force()
{
	for (int i = 0; i < npairs; i++) {
		delete pair[i];
		pair[i] = NULL;
	}
	memory->sfree(pair);
}

/* ---------------------------------------------------------------------- */

void Force::init()
{
	cut_max_global = 0.0;
	cut_min_global = BIG;

	int n = particle->ntypes;
	if (type2pair == NULL) {
		memory->create(type2pair, n+1, n+1, "Force: type2pair");
		for (int i = 0; i <= n; i++)
		for (int j = 0; j <= n; j++) {
			type2pair[i][j] = -1;
		}
	}

	for (int i = 0; i < npairs; i++) {
		pair[i]->init();
		cut_max_global = MAX(cut_max_global, pair[i]->cut_max);
		cut_min_global = MIN(cut_min_global, pair[i]->cut_min);
	}

	for (int i = 1; i <= n; i++)
	for (int j = 1; j <= n; j++) {
		if (type2pair[i][j] == -1) {
			char str[128];
			sprintf(str, "Pair type %d-%d has not been defined yet", i, j);
			error->warning(FLERR, str);
		}
	}
}

/* ---------------------------------------------------------------------- */

void Force::setup()
{
	for (int i = 0; i < npairs; i++) {
		pair[i]->setup();
	}
}

/* ----------------------------------------------------------------------
   Clear force
------------------------------------------------------------------------- */

void Force::clear()
{
	int nall = particle->nlocal + particle->nghost;
	//

	//for (int i = 0; i < nall; i++){
	//	particle->de[i] = 0.0;
	//	particle->drho[i] = 0.0;
	//	for (int j = 0; j < 3; j++) {
	//		particle->f[i][j] = 0.0;
	//		if (particle->torque_flag) particle->torque[i][j] = 0.0;
	//	}
	//}
	//

	//for (int i = 0; i < npairs; i++) 
	//for (int j = 0; j < 6; j++) {
	//	pair[i]->virial[j] = 0.0;
	//}
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	gpuforce_clear << < int(nall + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, BLOCK_SIZE >> >(nall, particle->devForceX, particle->devForceY, particle->devForceZ);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

}

/* ----------------------------------------------------------------------
   Compute force
------------------------------------------------------------------------- */

void Force::compute(int eflag, int vflag)
{
	// scan all pairs
	for (int i = 0; i < npairs; i++) {
		timer->stamp();
		pair[i]->compute(eflag, vflag);
		timer->stamp(TIME_PAIR1 + i);
	}
}

/* ----------------------------------------------------------------------
   Create pair style
------------------------------------------------------------------------- */

void Force::create_pair(int narg, char **arg)
{
	// create top level pair class
	if(npairs == maxpairs) {
		maxpairs += DELTA;
		pair = (Pair **) memory->srealloc(pair, maxpairs*sizeof(Pair *), "Force: pair");
	}
    
	if (0) return;

	// create second level pair class
#define PAIR_CLASS
#define PairStyle(key,Class) \
    else if (strcmp(arg[0],#key) == 0) \
		pair[npairs] = new Class(ps);
#include "style_pair.h"
#undef PairStyle
#undef PAIR_CLASS

	else error->all(FLERR, "Illegal pair style");

	pair[npairs]->set_style(narg,arg);            // set pair_style parameters in Pair Style class
	
	npairs++;
}

/* ----------------------------------------------------------------------
   Compute bounds implied by numeric str with a possible wildcard asterik
   1 = lower bound, nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = 1 to nmax,
     (3) i* = i to nmax, (4) *j = 1 to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void Force::bounds(char *str, int nmax, int &nlo, int &nhi)
{
	char *ptr = strchr(str,'*');

	if (ptr == NULL) {
		nlo = nhi = atoi(str);
	} else if (strlen(str) == 1) {
		nlo = 1;
		nhi = nmax;
	} else if (ptr == str) {
		nlo = 1;
		nhi = atoi(ptr+1);
	} else if (strlen(ptr+1) == 0) {
		nlo = atoi(str);
		nhi = nmax;
	} else {
		nlo = atoi(str);
		nhi = atoi(ptr+1);
	}

	if (nlo < 1 || nhi > nmax) error->all(FLERR,"Numeric index is out of bounds");
}

/* ----------------------------------------------------------------------
   find the pair id 
------------------------------------------------------------------------- */

int Force::find_pair_style(const char *pair_style)
{
	int pid;

	pid = -1;
	for (int i = 0; i < npairs; i++) {
		if (strcmp(pair_style, force->pair[i]->style) == 0) {	
			pid = i;       
		}
	}

	return pid;
}
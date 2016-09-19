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
#include "pointers.h"
#include "create_particle.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "parallel.h"
//#include "cuda_particle.h"
#include "cuda_sphparticle.h"
#include "particle.h"
#include "particle_type.h"
#include "phy_const.h"
#include "random_park.h"
#include "style_particle.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "device_atomic_functions.h"
#include "parallel.h"

using namespace PDPS_NS;
using namespace PhyConst;

#define DELTA 10000
#define EPSILON 1.0e-6

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

using namespace PDPS_NS;
using namespace std;

double CUDASPHParticle::LoadFactor = 0.5;
CUDASPHParticle::CUDASPHParticle(PDPS *ps) : Pointers(ps)
{
	HashTableSize = 0;
	//cudaAvec = NULL;
	devMapArray = NULL;
	devHashKey = NULL;
	devHashVal = NULL;
	devTry = NULL;
}

/* ---------------------------------------------------------------------- */

CUDASPHParticle::~CUDASPHParticle()
{
	if (devMapArray) cudaEngine->Free(devMapArray);
	if (devHashKey) cudaEngine->Free(devHashKey);
	if (devHashVal) cudaEngine->Free(devHashVal);
}

// atom map using map array
__global__ void gpuSetMap(
	int* __restrict Tag,
	int* __restrict MapArray,
	const int nAll
	)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < nAll) MapArray[Tag[i]] = i;
}

// atom map using hash table
__global__ void gpuHashMap(
	int*  __restrict Tag,
	unsigned int* __restrict HashKey,
	unsigned int* __restrict HashVal,
	const int N,
	const int nAll
	)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < nAll)
	{
		unsigned int p, now;
		unsigned int t = Tag[i];
		unsigned int u = t;
		unsigned int v = 0x4D6C3671;
		do {
			__TEA_CORE<8>(u, v);
			p = u % N;
			now = atomicCAS(HashKey + p, 0, t);
		} while (now && now != t);
		if (now == 0)
		{
			HashKey[p] = t;
			HashVal[p] = i;
		}
	}
}

void CUDASPHParticle::SetMap()
{
	int map_style = particle->map_style;
	int map_tag_max = particle->map_tag_max;
	int nlocal = particle->nlocal;
	int nghost = particle->nghost;
	int *tag = particle->tag;
	if (particle->map_style == 1){
		int max = 0;
		for (int i = 0; i < nlocal; i++) max = MAX(max, tag[i]);
		MPI_Allreduce(&max, &map_tag_max, 1, MPI_INT, MPI_MAX, mworld);

		if (cudaEngine->QueryMemSize(devMapArray) / sizeof(int) < map_tag_max + 1)
			cudaEngine->ReallocDevice("CUDAAtom::devMapArray", devMapArray, map_tag_max + 1, false, false);
		cudaMemsetAsync(devMapArray, 0XFF, (map_tag_max + 1) * sizeof(int), cudaEngine->StreamPool[0]);
		int threadsPerBlock = cudaEngine->QueryKernelBlockSize((void*)gpuSetMap);
	//	gpuSetMap << < n_block(nlocal + nghost, threadsPerBlock), threadsPerBlock, 0, cudaEngine->StreamPool[0] >> > (
	//		atom->cudaAvec->devTag,
	//		devMapArray,
	//		nlocal + nghost);
	}
	else
	{
		//	to be filled

	}

}


void CUDASPHParticle::data_atoms(int n, char *buf){

	//	to be filled

}
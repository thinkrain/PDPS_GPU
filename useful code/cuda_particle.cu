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
#include "cuda_particle.h"
#include "particle.h"
#include "particle_type.h"
#include "phy_const.h"
#include "random_park.h"
#include "style_particle.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "parallel.h"
#include "device_atomic_functions.h"

using namespace PDPS_NS;
using namespace PhyConst;

#define DELTA 10000
#define EPSILON 1.0e-6

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

using namespace PDPS_NS;
using namespace std;

double CUDAParticle::LoadFactor = 0.5;
CUDAParticle::CUDAParticle(PDPS *ps) : Pointers(ps)
{

	HashTableSize = 0;
	//cudaAvec = NULL;
	devMapArray = NULL;
	devHashKey = NULL;
	devHashVal = NULL;
	devTry = NULL;

	PermuKey = NULL;
	PermuTableA = NULL;
	PermuTableT = NULL;

	devTag = NULL;
	devType = NULL;
	devMask = NULL;
	devMass = NULL;
	devCoordX = NULL;
	devCoordY = NULL;
	devCoordZ = NULL;
	devForceX = NULL;
	devForceY = NULL;
	devForceZ = NULL;
	devVeloX = NULL;
	devVeloY = NULL;
	devVeloZ = NULL;
	devVirialXX = NULL;
	devVirialYY = NULL;
	devVirialZZ = NULL;
	devVirialXY = NULL;
	devVirialXZ = NULL;
	devVirialYZ = NULL;
	devImageX = NULL;
	devImageY = NULL;
	devImageZ = NULL;
	devEBond = NULL;
	devEAngle = NULL;
	devEDihed = NULL;
	devEImprop = NULL;
	devEPair = NULL;
	devRCoordX = NULL;
	devRCoordY = NULL;
	devRCoordZ = NULL;
	devRVeloX = NULL;
	devRVeloY = NULL;
	devRVeloZ = NULL;
	devMCoord = NULL;
	devMVelo = NULL;

	devHostCoord = NULL;
	devHostVelo = NULL;
	devHostForce = NULL;
	devHostImage = NULL;
	devHostTag = NULL;
	devHostType = NULL;
	devHostMask = NULL;
	devHostMassType = NULL;
	devHostSpecial = NULL;
	devHostMolecule = NULL;

	ptrHostCoord = NULL;
	ptrHostVelo = NULL;
	ptrHostForce = NULL;
	ptrHostImage = NULL;
	ptrHostTag = NULL;
	ptrHostType = NULL;
	ptrHostMask = NULL;
	ptrHostMassType = NULL;
	ptrHostSpecial = NULL;
	ptrHostMolecule = NULL;

	devHostNBond = devHostBondAtom = devHostBondType = NULL;
	ptrHostNBond = ptrHostBondAtom = ptrHostBondType = NULL;
	devHostNExcl = NULL;

	devNBond = NULL;
	devBonds = NULL;
	devBondsMapped = NULL;
	devNAngle = NULL;
	devAngles = NULL;
	devAngleType = NULL;
	devNDihed = NULL;
	devDiheds = NULL;
	devDihedType = NULL;
	devNImprop = NULL;
	devImprops = NULL;
	devImpropType = NULL;
	devNExcl = NULL;
	devExclTable = NULL;
}

/* ---------------------------------------------------------------------- */

CUDAParticle::~CUDAParticle()
{
	if (devMapArray) cudaEngine->Free(devMapArray);
	if (devHashKey) cudaEngine->Free(devHashKey);
	if (devHashVal) cudaEngine->Free(devHashVal);

	if (PermuKey)         cudaEngine->Free(PermuKey);
	if (PermuTableA)      cudaEngine->Free(PermuTableA);
	if (PermuTableT)      cudaEngine->Free(PermuTableT);

	if (devTag)           cudaEngine->Free(devTag);
	if (devType)          cudaEngine->Free(devType);
	if (devMask)          cudaEngine->Free(devMask);
	if (devMass)          cudaEngine->Free(devMass);
	if (devCoordX)        cudaEngine->Free(devCoordX);
	if (devCoordY)        cudaEngine->Free(devCoordY);
	if (devCoordZ)        cudaEngine->Free(devCoordZ);
	if (devForceX)        cudaEngine->Free(devForceX);
	if (devForceY)        cudaEngine->Free(devForceY);
	if (devForceZ)        cudaEngine->Free(devForceZ);
	if (devVeloX)         cudaEngine->Free(devVeloX);
	if (devVeloY)         cudaEngine->Free(devVeloY);
	if (devVeloZ)         cudaEngine->Free(devVeloZ);
	if (devVirialXX)      cudaEngine->Free(devVirialXX);
	if (devVirialYY)      cudaEngine->Free(devVirialYY);
	if (devVirialZZ)      cudaEngine->Free(devVirialZZ);
	if (devVirialXY)      cudaEngine->Free(devVirialXY);
	if (devVirialYZ)      cudaEngine->Free(devVirialYZ);
	if (devVirialXZ)      cudaEngine->Free(devVirialXZ);
	if (devImageX)        cudaEngine->Free(devImageX);
	if (devImageY)        cudaEngine->Free(devImageY);
	if (devImageZ)        cudaEngine->Free(devImageZ);
	if (devEBond)         cudaEngine->Free(devEBond);
	if (devEAngle)        cudaEngine->Free(devEAngle);
	if (devEDihed)        cudaEngine->Free(devEDihed);
	if (devEImprop)       cudaEngine->Free(devEImprop);
	if (devEPair)       cudaEngine->Free(devEPair);
	if (devRCoordX)       cudaEngine->Free(devRCoordX);
	if (devRCoordY)       cudaEngine->Free(devRCoordY);
	if (devRCoordZ)       cudaEngine->Free(devRCoordZ);
	if (devRVeloX)        cudaEngine->Free(devRVeloX);
	if (devRVeloY)        cudaEngine->Free(devRVeloY);
	if (devRVeloZ)        cudaEngine->Free(devRVeloZ);
	if (devMCoord)        cudaEngine->Free(devMCoord);
	if (devMVelo)        cudaEngine->Free(devMVelo);

	if (devNBond)       cudaEngine->Free(devNBond);
	if (devBonds)       cudaEngine->Free(devBonds);
	if (devBondsMapped) cudaEngine->Free(devBondsMapped);
	if (devNAngle)      cudaEngine->Free(devNAngle);
	if (devAngles)      cudaEngine->Free(devAngles);
	if (devAngleType)   cudaEngine->Free(devAngleType);
	if (devNDihed)      cudaEngine->Free(devNDihed);
	if (devDiheds)      cudaEngine->Free(devDiheds);
	if (devDihedType)   cudaEngine->Free(devDihedType);
	if (devNImprop)     cudaEngine->Free(devNImprop);
	if (devImprops)     cudaEngine->Free(devImprops);
	if (devImpropType) cudaEngine->Free(devImpropType);
	if (devNExcl)       cudaEngine->Free(devNExcl);
	if (devExclTable) cudaEngine->Free(devExclTable);
}

// particle map using map array
__global__ void gpuSetMap(
	int* __restrict Tag,
	int* __restrict MapArray,
	const int nAll
	)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < nAll) MapArray[Tag[i]] = i;
}

// particle map using hash table
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

void CUDAParticle::SetMap()
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
			cudaEngine->ReallocDevice("CUDAparticle::devMapArray", devMapArray, map_tag_max + 1, false, false);
		cudaMemsetAsync(devMapArray, 0XFF, (map_tag_max + 1) * sizeof(int), cudaEngine->StreamPool[0]);
		int threadsPerBlock = cudaEngine->QueryKernelBlockSize((void*)gpuSetMap);
	//	gpuSetMap << < n_block(nlocal + nghost, threadsPerBlock), threadsPerBlock, 0, cudaEngine->StreamPool[0] >> > (
	//		particle->cudaAvec->devTag,
	//		devMapArray,
	//		nlocal + nghost);
	}
	else
	{
		//	to be filled

	}

}


void CUDAParticle::data_particles(int n, char *buf){

	//	to be filled

}

void CUDAParticle::Grow(int nmaxOld, int nmaxNew)
{
#ifdef _RAINN_LOG_L3
	fprintf(stderr, "[CUDA] particle vector grew from %d to %d\n", nmaxOld, nmaxNew);
#endif

	cudaDeviceSynchronize();

	// per-particle array
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::PermuKey", PermuKey, nmaxNew, false, false);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::PermuTableA", PermuTableA, nmaxNew, false, false);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::PermuTableT", PermuTableT, nmaxNew, false, false);

	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devTag", devTag, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devType", devType, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devMask", devMask, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devMass", devMass, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devCoord", devCoordX, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devCoord", devCoordY, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devCoord", devCoordZ, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devForce", devForceX, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devForce", devForceY, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devForce", devForceZ, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devVelo", devVeloX, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devVelo", devVeloY, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devVelo", devVeloZ, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devVirial", devVirialXX, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devVirial", devVirialYY, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devVirial", devVirialZZ, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devVirial", devVirialXY, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devVirial", devVirialXZ, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devVirial", devVirialYZ, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devImage", devImageX, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devImage", devImageY, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devImage", devImageZ, nmaxNew);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devEBond", devEBond, nmaxNew, false, true);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devEAngle", devEAngle, nmaxNew, false, true);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devEDihed", devEDihed, nmaxNew, false, true);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devEImprop", devEImprop, nmaxNew, false, true);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devEPair", devEPair, nmaxNew, false, true);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devMCoord", devMCoord, nmaxNew, false, false);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devMVelo", devMVelo, nmaxNew, false, false);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devRCoord", devRCoordX, nmaxNew, false, false);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devRCoord", devRCoordY, nmaxNew, false, false);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devRCoord", devRCoordZ, nmaxNew, false, false);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devRVelo", devRVeloX, nmaxNew, false, false);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devRVelo", devRVeloY, nmaxNew, false, false);
	cudaEngine->ReallocDevice("CUDAparticleVecDPD::devRVelo", devRVeloZ, nmaxNew, false, false);

}

// pin host arries that hold primitive variables
void CUDAParticle::PinHostArray()
{
	cudaDeviceSynchronize();
	if (particle->nmax)
	{
		ptrHostCoord = cudaEngine->PinHost(devHostCoord, &(particle->x[0][0]), 3 * particle->nmax * sizeof(double));
		ptrHostForce = cudaEngine->PinHost(devHostForce, &(particle->f[0][0]), 3 * particle->nmax * sizeof(double));
		ptrHostVelo = cudaEngine->PinHost(devHostVelo, &(particle->v[0][0]), 3 * particle->nmax * sizeof(double));
		ptrHostType = cudaEngine->PinHost(devHostType, particle->type, particle->nmax * sizeof(int));
		ptrHostTag = cudaEngine->PinHost(devHostTag, particle->tag, particle->nmax * sizeof(int));
		ptrHostMask = cudaEngine->PinHost(devHostMask, particle->mask, particle->nmax * sizeof(int));
	}
}

// release host arries that hold primitive variables
void CUDAParticle::UnpinHostArray()
{
	cudaDeviceSynchronize();
	if (ptrHostCoord)    cudaEngine->UnpinHost(&(particle->x[0][0]), ptrHostCoord, devHostCoord);
	if (ptrHostForce)    cudaEngine->UnpinHost(&(particle->f[0][0]), ptrHostForce, devHostForce);
	if (ptrHostVelo)     cudaEngine->UnpinHost(&(particle->v[0][0]), ptrHostVelo, devHostVelo);
	if (ptrHostType)    cudaEngine->UnpinHost(particle->type, ptrHostType, devHostType);
	if (ptrHostTag)    cudaEngine->UnpinHost(particle->tag, ptrHostTag, devHostTag);
	if (ptrHostMask)    cudaEngine->UnpinHost(particle->mask, ptrHostMask, devHostMask);
	if (ptrHostMassType) cudaEngine->UnpinHost(particle->mass, ptrHostMassType, devHostMassType);

}
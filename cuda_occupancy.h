/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_CUDA_OCCUPANCY_H
#define PS_CUDA_OCCUPANCY_H

#include "pdps_cuda.h"
#include "pointers.h"

namespace PDPS_NS {

	struct OccuRecord
	{
		OccuRecord() { threadsPerBlock = 0, activeBlockPerSM = 0, occupancy = 0; }
		OccuRecord(int t, int b, float o) : threadsPerBlock(t), activeBlockPerSM(b), occupancy(o) {}
		int threadsPerBlock;
		int activeBlockPerSM;
		float occupancy;
	};

	struct GPUSPECS
	{
		int limitThreadsPerWarp;
		int limitWarpsPerMultiprocessor;
		int limitThreadsPerMultiprocessor;
		int limitBlocksPerMultiprocessor;
		int limitTotalSharedMemory;
		int limitTotalRegisters;
		int limitRegsPerThread;
		int maxThreadsPerBlock;
		int warpRegAllocGranularities;
		int myAllocationSize;
		int mySharedMemAllocationSize;
		int myWarpAllocationGranularity;
		int myAllocationGranularity;

		void read(string specs)
		{
			char strRegAllocGranu[32];
			sscanf(specs.c_str(),
				"%*f %*s %d	%d %d %d %d	%d %d %s %d	%d %d %d %*d %*d %*d %d",
				&limitThreadsPerWarp,
				&limitWarpsPerMultiprocessor,
				&limitThreadsPerMultiprocessor,
				&limitBlocksPerMultiprocessor,
				&limitTotalSharedMemory,
				&limitTotalRegisters,
				&myAllocationSize,
				strRegAllocGranu,
				&limitRegsPerThread,
				&mySharedMemAllocationSize,
				&myWarpAllocationGranularity,
				&maxThreadsPerBlock,
				&warpRegAllocGranularities
				);
			if (!strcmp(strRegAllocGranu, "block")) myAllocationGranularity = 1;
			else myAllocationGranularity = 0;
		}
	};

	class CUDAOccupancy : public Pointers
	{
	public:
		CUDAOccupancy(PDPS *);
		const vector<OccuRecord>& Occupancy(int device_id, void *kernel, size_t dynamic_shmem, enum cudaFuncCache shmemConfig);
		int2 ActiveBlock(int device_id, void *kernel, size_t dynamic_shmem, int threadsPerBlock, enum cudaFuncCache shmemConfig);
		int2 LeftPeak(int device_id, void *kernel, size_t dynamic_shmem, enum cudaFuncCache shmemConfig);
		int2 RightPeak(int device_id, void *kernel, size_t dynamic_shmem, enum cudaFuncCache shmemConfig);
		void ResolveSpecs();
	private:
		vector<string>   GPUSpecsTXT;
		vector<GPUSPECS> GPUSpecs;
		map<void *, vector<OccuRecord> > OccuRecTable;
	};

}

#endif
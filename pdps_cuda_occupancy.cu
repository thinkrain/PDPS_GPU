#include "pdps_cuda_engine.h"
#include "pdps_cuda_occupancy.h"
#include "pdps_cuda_util.h"

using namespace PDPS_NS;
using namespace std;

// initialize the GPU specification database
CUDAOccupancy::CUDAOccupancy(PDPS *ps) : Pointers(ps)
{
	GPUSpecsTXT.push_back( "1.0    sm_10    32    24     768    8    16384    8192    256    block    124    512    2     512    16384        0        0      0" );
	GPUSpecsTXT.push_back( "1.1    sm_11    32    24     768    8    16384    8192    256    block    124    512    2     512    16384        0        0      0" );
	GPUSpecsTXT.push_back( "1.2    sm_12    32    32    1024    8    16384    16384   512    block    124    512    2     512    16384        0        0      0" );
	GPUSpecsTXT.push_back( "1.3    sm_13    32    32    1024    8    16384    16384   512    block    124    512    2     512    16384        0        0      0" );
	GPUSpecsTXT.push_back( "2.0    sm_20    32    48    1536    8    49152    32768   64      warp     63    128    2    1024    49152    16384        0     64" );
	GPUSpecsTXT.push_back( "2.1    sm_21    32    48    1536    8    49152    32768   64      warp     63    128    2    1024    49152    16384        0     64" );
	GPUSpecsTXT.push_back( "3.0    sm_30    32    64    2048   16    49152    65536   256     warp     63    256    4    1024    49152    16384    32768    256" );
	GPUSpecsTXT.push_back( "3.5    sm_35    32    64    2048   16    49152    65536   256     warp    255    256    4    1024    49152    16384    32768    256" );
}

// parse the specification strings
void CUDAOccupancy::ResolveSpecs()
{
	if ( GPUSpecs.size() == cudaEngine->DevicePool.size() ) return;

	GPUSpecs.resize( cudaEngine->DevicePool.size() );
	for( int did = 0; did < GPUSpecs.size() ; did++ )
	{
		bool find = false;
		char sm_version[32];
		sprintf( sm_version, "%d.%d", cudaEngine->DeviceProperty[did].major, cudaEngine->DeviceProperty[did].minor );
		string smVersion = sm_version ;
		for( int i = 0 ; i < GPUSpecsTXT.size() ; i++ )
		{
			// matching specs for current architecture
			if ( GPUSpecsTXT[i].substr(0,3) == smVersion )
			{
				GPUSpecs[did].read(GPUSpecsTXT[i]);
				find = true;
				break;
			}
		}
		if ( !find )
		{
			fprintf(stderr,"[CDEV] cannot find matching specs for the current GPU #%d %s.\n", did, sm_version);
		}
	}
}

// calculate the occupancy of a given kernel at different block sizes
const vector<OccuRecord>& CUDAOccupancy::Occupancy( int device_id, void *kernel, size_t dynamic_shmem, enum cudaFuncCache shmemConfig )
{
	map<void *, vector<OccuRecord> >::iterator pRec = OccuRecTable.find( kernel );
	if ( pRec == OccuRecTable.end() )
	{
		vector<OccuRecord> Chart;

		ResolveSpecs();
		GPUSPECS &Specs = GPUSpecs[device_id];

		int kernelRegCount;
		int kernelSharedMemory;
		int configTotalSharedMemory;
		switch (shmemConfig)
		{
			case cudaFuncCachePreferL1    : configTotalSharedMemory = 16 * 1024 ; break;
			case cudaFuncCachePreferShared: configTotalSharedMemory = 48 * 1024 ; break;
			case cudaFuncCachePreferEqual : configTotalSharedMemory = 32 * 1024 ; break;
			default:                        configTotalSharedMemory = 48 * 1024 ; break;
		}

		// get kernel handle for querying attributes
		map<void *, cudaFuncAttributes >::iterator p = cudaEngine->KernelTable.find( kernel ) ;
		if ( p == cudaEngine->KernelTable.end() )
		{
			cudaFuncAttributes attr;
			cudaFuncGetAttributes( &attr, kernel );
			cudaEngine->KernelTable[ kernel ] = attr ;
			p = cudaEngine->KernelTable.find( kernel );
		}
		kernelRegCount = p->second.numRegs ;
		kernelSharedMemory = p->second.sharedSizeBytes + dynamic_shmem ;

		// calculate occupancy array with varying block size, fixed: reg usage & sharemem size
		Chart.reserve( Specs.maxThreadsPerBlock / Specs.limitThreadsPerWarp );
		for( int MyThreadCount = Specs.limitThreadsPerWarp ; MyThreadCount <= Specs.maxThreadsPerBlock ; MyThreadCount += Specs.limitThreadsPerWarp )
		{
			// CUDA Occupancy Calculator, B38-B40
			int MyWarpsPerBlock, MyRegsPerBlock, MySharedMemPerBlock ;
			MyWarpsPerBlock = ceiling( MyThreadCount, Specs.limitThreadsPerWarp ) / Specs.limitThreadsPerWarp ;
			if ( Specs.myAllocationGranularity == 1 )
				MyRegsPerBlock = ceiling( ceiling( MyWarpsPerBlock, Specs.myWarpAllocationGranularity ) * kernelRegCount * Specs.limitThreadsPerWarp, Specs.myAllocationSize );
			else
				MyRegsPerBlock = ceiling( kernelRegCount * Specs.limitThreadsPerWarp, Specs.myAllocationSize ) * ceiling( MyWarpsPerBlock, Specs.myWarpAllocationGranularity );
			MySharedMemPerBlock = kernelSharedMemory ;

			// CUDA Occupancy Calculator, D38-D40
			int limitBlocksDueToWarps, limitBlocksDueToRegs, limitBlocksDueToSMem;
			limitBlocksDueToWarps = min( Specs.limitBlocksPerMultiprocessor, Specs.limitWarpsPerMultiprocessor / MyWarpsPerBlock ) ;
			if ( kernelRegCount > Specs.limitRegsPerThread ) limitBlocksDueToRegs = 0;
			else
			{
				if ( kernelRegCount > 0 ) limitBlocksDueToRegs = Specs.limitTotalRegisters / MyRegsPerBlock ;
				else limitBlocksDueToRegs = Specs.limitBlocksPerMultiprocessor ;
			}
			if ( MySharedMemPerBlock > 0 ) limitBlocksDueToSMem = configTotalSharedMemory / MySharedMemPerBlock ;
			else limitBlocksDueToSMem = Specs.limitBlocksPerMultiprocessor ;

			// Calculate occupancy
			int ActiveBlocks  = min( min( limitBlocksDueToWarps, limitBlocksDueToRegs ), limitBlocksDueToSMem ) ;
			int ActiveWarps   = ActiveBlocks * MyWarpsPerBlock ;
			int ActiveThreads = ActiveWarps  * Specs.limitThreadsPerWarp ;
			float Occupancy  = (double)ActiveWarps / Specs.limitWarpsPerMultiprocessor ;
			Chart.push_back( OccuRecord( MyThreadCount, ActiveBlocks, Occupancy ) );
		}

		OccuRecTable[ kernel ] = Chart;
		pRec = OccuRecTable.find( kernel );
	}
	return pRec->second;
}

// decide the number of active blocks per SM for a given launch config
// .x: active block per SM
// .y: active block on entire GPU
int2 CUDAOccupancy::ActiveBlock( int device_id, void *kernel, size_t dynamic_shmem, int threadsPerBlock, enum cudaFuncCache shmemConfig )
{
	ResolveSpecs();
	GPUSPECS &Specs = GPUSpecs[device_id];

	int kernelRegCount;
	int kernelSharedMemory;
	int configTotalSharedMemory;
	switch (shmemConfig)
	{
		case cudaFuncCachePreferL1    : configTotalSharedMemory = 16 * 1024 ; break;
		case cudaFuncCachePreferShared: configTotalSharedMemory = 48 * 1024 ; break;
		case cudaFuncCachePreferEqual : configTotalSharedMemory = 32 * 1024 ; break;
		default:                        configTotalSharedMemory = 48 * 1024 ; break;
	}

	// get kernel handle for querying attributes
	map<void *, cudaFuncAttributes >::iterator p = cudaEngine->KernelTable.find( kernel ) ;
	if ( p == cudaEngine->KernelTable.end() )
	{
		cudaFuncAttributes attr;
		cudaFuncGetAttributes( &attr, kernel );
		cudaEngine->KernelTable[ kernel ] = attr ;
		p = cudaEngine->KernelTable.find( kernel );
	}
	kernelRegCount = p->second.numRegs ;
	kernelSharedMemory = p->second.sharedSizeBytes + dynamic_shmem ;

	// CUDA Occupancy Calculator, B38-B40
	int MyWarpsPerBlock, MyRegsPerBlock, MySharedMemPerBlock ;
	MyWarpsPerBlock = ceiling( threadsPerBlock, Specs.limitThreadsPerWarp ) / Specs.limitThreadsPerWarp ;
	if ( Specs.myAllocationGranularity == 1 )
		MyRegsPerBlock = ceiling( ceiling( MyWarpsPerBlock, Specs.myWarpAllocationGranularity ) * kernelRegCount * Specs.limitThreadsPerWarp, Specs.myAllocationSize );
	else
		MyRegsPerBlock = ceiling( kernelRegCount * Specs.limitThreadsPerWarp, Specs.myAllocationSize ) * ceiling( MyWarpsPerBlock, Specs.myWarpAllocationGranularity );
	MySharedMemPerBlock = kernelSharedMemory ;

	// CUDA Occupancy Calculator, D38-D40
	int limitBlocksDueToWarps, limitBlocksDueToRegs, limitBlocksDueToSMem;
	limitBlocksDueToWarps = min( Specs.limitBlocksPerMultiprocessor, Specs.limitWarpsPerMultiprocessor / MyWarpsPerBlock ) ;
	if ( kernelRegCount > Specs.limitRegsPerThread ) limitBlocksDueToRegs = 0;
	else
	{
		if ( kernelRegCount > 0 ) limitBlocksDueToRegs = Specs.limitTotalRegisters / MyRegsPerBlock ;
		else limitBlocksDueToRegs = Specs.limitBlocksPerMultiprocessor ;
	}
	if ( MySharedMemPerBlock > 0 ) limitBlocksDueToSMem = configTotalSharedMemory / MySharedMemPerBlock ;
	else limitBlocksDueToSMem = Specs.limitBlocksPerMultiprocessor ;

	// Calculate occupancy
	int Active = min( min( limitBlocksDueToWarps, limitBlocksDueToRegs ), limitBlocksDueToSMem ) ;
	return make_int2( Active, Active * cudaEngine->DeviceProperty[0].multiProcessorCount );
}

// determine the kernel config with maximum occupancy and largest block size
int2 CUDAOccupancy::LeftPeak( int device_id, void *kernel, size_t dynamic_shmem, enum cudaFuncCache shmemConfig )
{
	const vector<OccuRecord> &Chart = Occupancy( device_id, kernel, dynamic_shmem, shmemConfig );

	float peakOccu = 0.;
	for( int i = 0 ; i < Chart.size() ; i++ )
		peakOccu = max( peakOccu, Chart[i].occupancy ) ;
	for( vector<OccuRecord>::const_iterator i = Chart.begin() ; i != Chart.end() ; i++ )
	{
		if ( fabs( i->occupancy - peakOccu ) < 1e-4 )
		{
			int2 ret;
			ret.x = i->activeBlockPerSM * cudaEngine->DeviceProperty[ device_id ].multiProcessorCount;
			ret.y = i->threadsPerBlock ;
			return ret;
		}
	}
	int2 ret;
	ret.x = Chart.back().activeBlockPerSM * cudaEngine->DeviceProperty[ device_id ].multiProcessorCount;
	ret.y = Chart.back().threadsPerBlock;
	return ret;
}

// determine the kernel config with maximum occupancy and smallest block size
int2 CUDAOccupancy::RightPeak( int device_id, void *kernel, size_t dynamic_shmem, enum cudaFuncCache shmemConfig )
{
	const vector<OccuRecord> &Chart = Occupancy( device_id, kernel, dynamic_shmem, shmemConfig );

	float peakOccu = 0.;
	for( int i = 0 ; i < Chart.size() ; i++ )
		peakOccu = max( peakOccu, Chart[i].occupancy ) ;
	for( vector<OccuRecord>::const_reverse_iterator i = Chart.rbegin() ; i != Chart.rend() ; i++ )
	{
		if ( fabs( i->occupancy - peakOccu ) < 1e-4 )
		{
			int2 ret;
			ret.x = i->activeBlockPerSM * cudaEngine->DeviceProperty[ device_id ].multiProcessorCount ;
			ret.y = i->threadsPerBlock ;
			return ret;
		}
	}
	int2 ret;
	ret.x = Chart.back().activeBlockPerSM * cudaEngine->DeviceProperty[ device_id ].multiProcessorCount ;
	ret.y = Chart.back().threadsPerBlock ;
	return ret;
}



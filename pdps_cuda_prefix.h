#ifndef _PS_CUDA_PREFIX
#define _PS_CUDA_PREFIX

// specialized prefix sum kernel
// 1G bit prefix sum by 3 passes

#include "pdps_cuda.h"
#include "pdps_cuda_engine.h"

using namespace std;

namespace LAMMPS_NS
{

// digest arbitrary array into at most 1024 sectors
template<class TYPE> __global__ void gpuBitPrefixDown( uint *in, TYPE *out, TYPE *next, const int jobsize, const int njob, const int n )
{
	__shared__ int previous_sum;
	__shared__ int buffer[32];

	for(int j = blockIdx.x; j < njob; j += gridDim.x)
	{
		if ( threadIdx.x == 0 ) previous_sum = 0;
		__syncthreads();

		for(int i = threadIdx.x ; i < jobsize ; i += blockDim.x)
		{
			// ld
			int p = i + j * jobsize;
			int v = ( p < n ) ? ( __popc(in[p]) ) : 0;

			// intra-warp prefix
			int prefix_warp = __warp_prefix_excl( v );
			if ( __laneid() == 31 ) buffer[ threadIdx.x / 32 ] = prefix_warp + v ; // + x : inclusive result
			__syncthreads();

			// inter-warp prefix
			if ( threadIdx.x < 32 ) buffer[ threadIdx.x ] = __warp_prefix_excl( buffer[threadIdx.x] );
			__syncthreads();

			// intra-warp shift & st
			int prefix_block = prefix_warp + buffer[ threadIdx.x / 32 ];
			if ( p < n ) out[p] = prefix_block + previous_sum ;
			__syncthreads();
			if ( threadIdx.x == blockDim.x - 1 ) previous_sum += prefix_block + v ;
		}
		__syncthreads();

		if ( threadIdx.x == 0) next[ j ] = previous_sum;
		__syncthreads();
	}
}

// single CTA 1024 element prefix sum
template<class TYPE> __global__ void gpuBitPrefixBlock( TYPE *val, const int n )
{
	__shared__ int buffer[32];

	// ld & intra-warp prefix
	int v = ( threadIdx.x < n ) ? ( val[threadIdx.x] ) : 0;
	int prefix_warp = __warp_prefix_excl( v );
	if ( __laneid() == 31 ) buffer[ threadIdx.x / 32 ] = prefix_warp + v ; // + x : inclusive result
	__syncthreads();

	// inter-warp prefix
	if ( threadIdx.x < 32 ) buffer[ threadIdx.x ] = __warp_prefix_excl( buffer[threadIdx.x] );
	__syncthreads();

	// intra-warp shift & st
	if ( threadIdx.x < n ) val[threadIdx.x] = prefix_warp + buffer[ threadIdx.x / 32 ];
	__syncthreads();
}

// broadcast the single CTA result back to original array
template<class TYPE> __global__ void gpuBitPrefixUp( TYPE *out, TYPE *next, const int work, const int n )
{
	for(int i = blockIdx.x * blockDim.x + threadIdx.x ; i < n ; i+= gridDim.x * blockDim.x )
	{
		out[i] += next[ i / work ];
	}
}

template<class TYPE> class CUDAPrefixPlan : protected Pointers
{
public:
	CUDAPrefixPlan(LAMMPS *lmp) : Pointers(lmp)
	{
		L1 = NULL;
		L2 = NULL;
		GridCfg1 = make_int2( 0, 0 );
		GridCfg2 = make_int2( 0, 0 );
		GridCfg3 = make_int2( 0, 0 );
		CreatePlan( 1024*1024 );
	}
	// copy constructor
	CUDAPrefixPlan( const CUDAPrefixPlan<TYPE> &another ) : Pointers(another.lmp) // copy constructor
	{
		L1 = NULL;
		L2 = NULL;
		GridCfg1 = make_int2( 0, 0 );
		GridCfg2 = make_int2( 0, 0 );
		GridCfg3 = make_int2( 0, 0 );
		CreatePlan( 1024*1024 );
	}
	// copy-assign operator
	CUDAPrefixPlan<TYPE>& operator = ( const CUDAPrefixPlan<TYPE> &another ) // copy assignment operator
	{
		lmp      = another.lmp;
		memory   = another.lmp->memory;
		error    = another.lmp->error;
		universe = another.lmp->universe;
		input    = another.lmp->input;
		atom     = another.lmp->atom;
		update   = another.lmp->update;
		neighbor = another.lmp->neighbor;
		comm     = another.lmp->comm;
		domain   = another.lmp->domain;
		force    = another.lmp->force;
		modify   = another.lmp->modify;
		group    = another.lmp->group;
		output   = another.lmp->output;
		timer    = another.lmp->timer;
		world    = another.lmp->world;
		infile   = another.lmp->infile;
		screen   = another.lmp->screen;
		logfile  = another.lmp->logfile;
		cudaEngine = another.lmp->cudaEngine;
		DestroyPlan();
		CreatePlan( another.N );
		return *this;
	}
	~CUDAPrefixPlan() { DestroyPlan(); }
	// perform the prefix sum on the specified stream
	int* Prefix( uint *Bits, int N, cudaStream_t stream = 0 )
	{
		if ( this->N < N )
		{
			cudaStreamSynchronize( stream );
			DestroyPlan();
			CreatePlan( N * 1.2 );
		}
		GridConfig();
		int nInt = ceiling(N,32) / 32;
		int job_size = ( nInt < GridCfg1.y * GridCfg2.y ) ? ( GridCfg1.y ) : ceiling( nInt / GridCfg2.y, GridCfg1.y ) ;
		int njob = n_block( nInt, job_size );
		gpuBitPrefixDown <<< GridCfg1.x, GridCfg1.y, 0, stream >>>( Bits, L1, L2, job_size, njob, nInt );
		gpuBitPrefixBlock<<<          1, GridCfg2.y, 0, stream >>>( L2, n_block( nInt, job_size ) );
		gpuBitPrefixUp   <<< GridCfg3.x, GridCfg3.y, 0, stream >>>( L1, L2, job_size, nInt );
		return L1;
	}
public:
	int N;
	TYPE *L1, *L2;
	int2 GridCfg1, GridCfg2, GridCfg3;

	void CreatePlan( int N ) // N: number of bits
	{
		if (!N) return;
		this->N		   = N;
		cudaEngine->ReallocDevice( "CUDAPrefixPlan::L1", L1, ceiling(N,32) / 32, false, false );
		cudaEngine->ReallocDevice( "CUDAPrefixPlan::L2", L2, 1024, false, false );
	}
	void DestroyPlan()
	{
		N = 0;
		if ( L1 ) cudaEngine->Free( L1 );
		if ( L2 ) cudaEngine->Free( L2 );
	}
	// determine optimal kernel configuration on the GPU
	void GridConfig()
	{
		if ( !GridCfg1.x )
		{
			GridCfg1 = cudaEngine->OccuCalc.RightPeak( 0, (void*)gpuBitPrefixDown <TYPE>, 0, cudaFuncCachePreferShared );
			GridCfg2 = cudaEngine->OccuCalc.RightPeak( 0, (void*)gpuBitPrefixBlock<TYPE>, 0, cudaFuncCachePreferShared );
			GridCfg3 = cudaEngine->OccuCalc.RightPeak( 0, (void*)gpuBitPrefixUp   <TYPE>, 0, cudaFuncCachePreferShared );
		}
	}
};

}

#endif


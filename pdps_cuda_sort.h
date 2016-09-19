#ifndef _PS_CUDA_SORT
#define _PS_CUDA_SORT

#include "pdps_cuda.h"
#include "pdps_cuda_engine.h"

using namespace std;

namespace PDPS_NS
{

// first stage of a radix sort pass: extract radices and perform prefix sum
template<typename TYPE> __global__ void gpuRadixPrefixSum( TYPE *vals, const int n, const int padding )
{
	#ifdef _RAINN_DEBUG
	assert( warpSize == 32 );
	#endif

	__shared__ TYPE previous_sum;
	__shared__ TYPE buffer[32];
	if ( threadIdx.x == 0 ) previous_sum = 0;
	__syncthreads();

	int len = ( n + blockDim.x - 1 ) & ( ~( blockDim.x - 1 ) );
	TYPE *val = vals + blockIdx.x * padding;
	for(int i = threadIdx.x ; i < len ; i+= blockDim.x )
	{
		// ld
		TYPE in;
		if ( i < n ) in = val[i];
		else in = 0;

		// intra-warp prefix
		TYPE sum_warp = __warp_prefix_excl( in );
		if ( __laneid() == 31 ) buffer[ threadIdx.x / 32 ] = sum_warp + in ; // + x : inclusive result
		__syncthreads();

		// inter-warp prefix
		if ( threadIdx.x < 32 ) buffer[ threadIdx.x ] = __warp_prefix_excl( buffer[threadIdx.x] );
		__syncthreads();

		// intra-warp shift & st
		TYPE sum_block = sum_warp + buffer[ threadIdx.x / 32 ];
		if ( i < n ) val[i] = sum_block + previous_sum ;
		__syncthreads();
		if ( threadIdx.x == blockDim.x - 1 ) previous_sum += sum_block + in ;
	}
	__syncthreads();
	if ( threadIdx.x == 0) val[padding-1] = previous_sum;
}

// second stage of a radix sort pass: build histogram of each radix
template<int RADIX, typename TYPE>
void __global__ gpuRadixHistogram(
	TYPE* __restrict key,
	int*  __restrict tprev,
	int*  __restrict wprev,
	const int len,
	const int cur_bit,
	const int prev_padding
)
{
	__shared__ int *wLocal[RADIX];
	if ( threadIdx.x < RADIX ) wLocal[ threadIdx.x ] = wprev + threadIdx.x * prev_padding ;
	__syncthreads();

	int laneID = __laneid();
	for( int i = threadIdx.x + blockIdx.x * blockDim.x ; i < len ; i += blockDim.x * gridDim.x )
	{
		int radix = ( key[i] >> cur_bit ) & (RADIX-1);
		#pragma unroll
		for( int r = 0 ; r < RADIX ; r++ )
		{
			uint before = __ballot( radix == r );
			if ( radix  == r ) tprev[i] = __popc( before << ( WARPSZ - laneID ) );
			if ( laneID == 0 ) wLocal[ r ][ i / WARPSZ ] = __popc(before) ;
	   	}
	}
}

// third stage of a radix sort pass: scatter the key/values pairs based on the result of prefix sum and histogram
template<int RADIX,
		 typename KeyType,
		 typename ValType>
void __global__ gpuRadixPermute(
	KeyType* __restrict key,
	ValType* __restrict val,
	KeyType* __restrict key_out,
	ValType* __restrict val_out,
		int* __restrict tprev,
		int* __restrict wprev,
	const int len,
	const int cur_bit,
	const int prev_padding
)
{
	// get_num_of_total_ones
	__shared__ int totals[RADIX];
	if ( threadIdx.x == 0 )
	{
		totals[0] = 0;
		for( int p = 1 ; p < RADIX ; p++ ) totals[p] = totals[p-1] + wprev[ p * prev_padding - 1 ];
	}
	__syncthreads();

	for( int i = threadIdx.x + blockIdx.x * blockDim.x ; i < len ; i += blockDim.x * gridDim.x )
	{
		int radix = ( key[i] >> cur_bit ) & (RADIX-1);
		int new_i = totals[radix] + tprev[i] + wprev[ radix * prev_padding+ (i/WARPSZ) ] ;
		key_out[ new_i ] = key[i];
		val_out[ new_i ] = val[i];
	}
}

inline double LOG2(double x){ return log(x)/log(2.0); }

// this is the template class for radix sort on GPU supporting arbitrary data types
template<int RADIX, class KeyType, class ValType> class CUDASortPlan : Pointers
{
public:
	CUDASortPlan(LAMMPS *lmp) : Pointers(lmp)
	{
		CreatePlan( 1024 );
		GridCfg1 = make_int2(0,0);
		GridCfg2 = make_int2(0,0);
		GridCfg3 = make_int2(0,0);
	}
	~CUDASortPlan() { DestroyPlan(); }
	// sort with nBit per pass
	void Sort( KeyType *Key, ValType *Val, int N, int nBit, cudaStream_t stream = 0 )
	{
		if ( this->N < N || this->N > N * 2 )
		{
			cudaStreamSynchronize( stream );
			DestroyPlan();
			CreatePlan( N * 1.2 );
		}
		GridConfig();
		KeyType *KeyIn = Key;
		ValType *ValIn = Val;
		for(int bit = 0; bit < nBit; bit += LOG2(RADIX) )
		{
			gpuRadixHistogram<RADIX, KeyType>       <<< GridCfg1.x, GridCfg1.y, 0, stream>>>( KeyIn, ThreadPrev, WarpPrev, N, bit, WarpPadding );
			gpuRadixPrefixSum                       <<<      RADIX, GridCfg2.y, 0, stream>>>( WarpPrev, ( N + WARPSZ - 1 ) / WARPSZ, WarpPadding);
			gpuRadixPermute<RADIX, KeyType, ValType><<< GridCfg3.x, GridCfg3.y, 0, stream>>>( KeyIn, ValIn, KeyOut, ValOut, ThreadPrev, WarpPrev, N, bit, WarpPadding );
			swap( KeyIn, KeyOut );
			swap( ValIn, ValOut );
		}
		if ( Key != KeyIn )
		{
			cudaMemcpyAsync( Key, KeyIn, sizeof(KeyType)*N, cudaMemcpyDefault, stream );
			cudaMemcpyAsync( Val, ValIn, sizeof(ValType)*N, cudaMemcpyDefault, stream );
			swap( KeyIn, KeyOut );
			swap( ValIn, ValOut );
		}
	}
protected:
	int N, WarpPadding;
	int *ThreadPrev, *WarpPrev; // radix counting buffer
	KeyType *KeyOut;
	ValType *ValOut; // ping-pong buffer
	int2 GridCfg1, GridCfg2, GridCfg3;

	// create buffers
	void CreatePlan( int N )
	{
		if (!N) return;
		this->N		   = N;
		this->WarpPadding = ( N + WARPSZ - 1 ) / WARPSZ + 1;
		this->KeyOut	  = cudaEngine->MallocDevice<KeyType>( "CUDASortPlan::KeyOut", N );
		this->ValOut	  = cudaEngine->MallocDevice<ValType>( "CUDASortPlan::ValOut", N );
		this->ThreadPrev  = cudaEngine->MallocDevice<int>( "CUDASortPlan::ThreadPrev", N );
		this->WarpPrev	  = cudaEngine->MallocDevice<int>( "CUDASortPlan::WarpPrev", WarpPadding * RADIX );
	}
	// free buffers
	void DestroyPlan()
	{
		N = WarpPadding = 0;
		if ( KeyOut )	  cudaEngine->Free( KeyOut );
		if ( ValOut )	  cudaEngine->Free( ValOut );
		if ( ThreadPrev ) cudaEngine->Free( ThreadPrev );
		if ( WarpPrev )   cudaEngine->Free( WarpPrev );
	}
	// determine optimal GPU kernel configuration
	void GridConfig()
	{
		if ( !GridCfg1.x )
		{
			GridCfg1 = cudaEngine->OccuCalc.RightPeak( 0, (void*)gpuRadixHistogram<RADIX,KeyType>,         0, cudaFuncCachePreferShared );
			GridCfg2 = cudaEngine->OccuCalc.RightPeak( 0, (void*)gpuRadixPrefixSum<int>,                0, cudaFuncCachePreferShared );
			GridCfg3 = cudaEngine->OccuCalc.RightPeak( 0, (void*)gpuRadixPermute<RADIX,KeyType,ValType>, 0, cudaFuncCachePreferShared );
		}
	}
};

}

#endif


#ifndef _PS_CUDA_UTIL
#define _PS_CUDA_UTIL

#include <cmath>
#include <assert.h>
#include "cuda.h"
#include "cuda_runtime.h"
#include "math_constants.h"

using namespace std;

#define EPSILON         1.0E-10
#define EPSILON_SQ      1.0E-20
#define _LN_2           6.9314718055994528623E-1
#define _LOG2_E         1.4426950408889634074
#define _1_OVER_SQ2     7.0710678118654757274E-1
#define _SQRT_2         1.4142135623730950488
#define _SQRT_3         1.7320508075688772935
#define _2_TO_MINUS_30  9.3132257461547851562E-10
#define _2_TO_MINUS_31  4.6566128730773925781E-10
#define _2_TO_MINUS_32  2.3283064365386962891E-10

#if(__CUDA_ARCH__<=350)
#define WARPSZ 32
#else
#error UNKNOWN ARCHITECTURE FOR DETERMINING WARP SIZE
#endif

inline int ceiling( int x, int inc )
{
	return ( ( x + inc - 1 ) / inc ) * inc ;
}

#define _TEA_K0      0xA341316C
#define _TEA_K1      0xC8013EA4
#define _TEA_K2      0xAD90777D
#define _TEA_K3      0x7E95761E
#define _TEA_DT      0x9E3779B9

template<int N> __inline__ __host__ __device__ void __TEA_CORE(unsigned int &v0, unsigned int &v1, unsigned int sum = 0)
{
	sum += _TEA_DT;
	v0 += ((v1 << 4) + _TEA_K0) ^ (v1 + sum) ^ ((v1 >> 5) + _TEA_K1);
	v1 += ((v0 << 4) + _TEA_K2) ^ (v0 + sum) ^ ((v0 >> 5) + _TEA_K3);
	__TEA_CORE< N - 1 >(v0, v1, sum);
}

template<> __inline__ __host__ __device__ void __TEA_CORE<1>(unsigned int &v0, unsigned int &v1, unsigned int sum)
{
	sum += _TEA_DT;
	v0 += ((v1 << 4) + _TEA_K0) ^ (v1 + sum) ^ ((v1 >> 5) + _TEA_K1);
	v1 += ((v0 << 4) + _TEA_K2) ^ (v0 + sum) ^ ((v0 >> 5) + _TEA_K3);
}

// calculate the number of blocks required to hold given number of threads
inline unsigned int n_block( unsigned int threads, unsigned int threadPerBlock )
{
	return ( threads + threadPerBlock - 1 ) / threadPerBlock;
}

#endif

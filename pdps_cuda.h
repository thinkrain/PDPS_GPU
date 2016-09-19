#ifndef PS_CUDA
#define PS_CUDA

#define DEBUG_TAG0 171484
#define DEBUG_TAG1 171485
#define DEBUG_TAG2 171486

#include <cmath>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>

#include "omp.h"
#include "mpi.h"
#include "cuda.h"
#include "pdps_cuda_util.h"
#include "cuda_runtime.h"

#include "pdps.h"
#include "pointers.h"
//#include "cuda_util.h"
#include "cuda_profiler_api.h"

using namespace std;

namespace LAMMPS_NS {

//#define _RAINN_LOG_MEM

#define MAX_BLOCK_SIZE 1024

#define CUDA_DEVICE_NUM_ZERO											0xAE000001
#define CUDA_DEVICE_LIST_EMPTY											0xAE000002
#define CUDA_ENGINE_ALLOCATION_FAIL										0xAE000003
#define CUDA_DEVICE_CANNOT_MAP_HOST_MEMORY								0xAE000004
#define CUDA_DEVICE_NOT_SUPPORT_UVA										0xAE000005
#define OMP_NUM_THREADS_TOO_SMALL										0xAE000006
#define CUDA_BOND_NUMBER_EXCEED_LIMIT									0xAE100001
#define CUDA_KERNEL_CACHE_PREF_NOT_SET									0xAE200001

typedef double	     r64;
typedef float	     r32;
typedef float4       r32v;
typedef unsigned long long u64;

#define CUDAVEC_TYPEMASSTAG  1
#define CUDAVEC_COORD        2
#define CUDAVEC_FORCE        4
#define CUDAVEC_VELO         8
#define CUDAVEC_IMAGE       16
#define CUDAVEC_EXCL        32
#define CUDAVEC_BOND        64
#define CUDAVEC_MOLE       128
// to be continued...

#define CUDAVEC_LOCAL        1
#define CUDAVEC_GHOST        2

#define CUDACOPY_C2G         0
#define CUDACOPY_G2C         1

// all: from init to finalize
// loop: from integrate::setup to last step of integrate::iterate
// core: from the 25% - 75% steps
#define CUDAPROF_NONE        0
#define CUDAPROF_ALL         1
#define CUDAPROF_LOOP        2
#define CUDAPROF_CORE        3

#define _FAST_DEBUG_	cudaDeviceSynchronize();\
						cout<<cudaGetErrorString(cudaGetLastError())<<' '<<__LINE__<<' '<<__FILE__<<endl;
#define _FAST_EXIT_		cudaDeviceSynchronize();\
						cout<<cudaGetErrorString(cudaGetLastError())<<' '<<__LINE__<<' '<<__FILE__<<endl;\
						cudaDeviceReset();\
						exit(0);
}

#endif

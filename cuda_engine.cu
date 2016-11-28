/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
//#include "stdint.h"

#include "error.h"
#include "cuda_engine.h"

using namespace PDPS_NS;
using namespace std;

CUDAEngine::CUDAEngine(PDPS *ps, int Device, string Profile) : Pointers(ps), OccuCalc(ps)
{
	Init(Device, Profile);
}

CUDAEngine::~CUDAEngine()
{
#ifdef _RAINN_LOG_L3
	cout << "Kernel Summary:" << endl;
	map<void *, cudaFuncAttributes >::iterator p;
	for (p = KernelTable.begin(); p != KernelTable.end(); p++)
	{
		fprintf(stderr, "Kernel    : %p \n", p->first);
		fprintf(stderr, "numRegs   : %d \n", p->second.numRegs);
		fprintf(stderr, "SharedMem : %.2f KB \n", p->second.sharedSizeBytes / 1024.0);
		fprintf(stderr, "ConstMem  : %.2f KB \n", p->second.constSizeBytes / 1024.0);
		fprintf(stderr, "LocalMem  : %.2f KB \n", p->second.localSizeBytes / 1024.0);
		fprintf(stderr, "BlockSize : %d \n", p->second.maxThreadsPerBlock);
	}
#endif

	Destroy();
}

// init the management layer
// fetch device, allocate streams
int CUDAEngine::Init(int Device, string Profile)
{
	if (Profile == "all") ProfilingMode = CUDAPROF_ALL;
	else if (Profile == "loop") ProfilingMode = CUDAPROF_LOOP;
	else if (Profile == "core") ProfilingMode = CUDAPROF_CORE;
	else ProfilingMode = CUDAPROF_NONE;
	// switch to the device
	DevicePool.push_back(Device);
	cudaSetDevice(Device);
	// start profile if in CUDAPROF_ALL mode
	if (ProfilingMode == CUDAPROF_ALL) StartProfiler();
	// check capabilities & logging
	cudaDeviceProp DeviceProp;
	cudaGetDeviceProperties(&DeviceProp, Device);
	DeviceProperty.push_back(DeviceProp);
	if (!DeviceProp.canMapHostMemory)
	{
		fprintf(stderr, "ERR 0X%X, Mapped memory is not supported by this GPU, quiting...", CUDA_DEVICE_CANNOT_MAP_HOST_MEMORY);
		exit(CUDA_DEVICE_CANNOT_MAP_HOST_MEMORY);
	}
	//if (!DeviceProp.unifiedAddressing)
	//{
	//	fprintf(stderr, "ERR 0X%X, Unified Virtual Addressing (UVA) is not supported by this GPU, quiting...", CUDA_DEVICE_NOT_SUPPORT_UVA);
	//	exit(CUDA_DEVICE_NOT_SUPPORT_UVA);
	//}
	// create stream
	for (int s = 0; s < 16; s++)
	{
		cudaStream_t stream;
		cudaStreamCreate(&stream);
		StreamPool.push_back(stream);
	}
	cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
	cudaDeviceSetLimit(cudaLimitMallocHeapSize, DeviceProp.totalGlobalMem * 0.2);
	cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 32 * 1024 * 1024);
	// display information
#ifdef _RAINN_LOG_L1
#ifdef _RAINN_LOG_L2
	fprintf(stderr, "[CUDA] \tDevice %d, Compute capability: %d.%d\n", Device, DeviceProp.major, DeviceProp.minor);
	fprintf(stderr, "[CUDA] \tRunning at %.2f GHz, L2 cache size = %d kB\n", DeviceProp.clockRate / 1048576.0, DeviceProp.l2CacheSize / 1024);
	fprintf(stderr, "[CUDA] \tGRAM size = %d MB bandwidth = %d bit running at %d MHz\n", DeviceProp.totalGlobalMem / 1024 / 1024, DeviceProp.memoryBusWidth, DeviceProp.memoryClockRate / 1024);
#endif
	for (int s = 0; s < 16; s++)
		fprintf(stderr, "[CUDA] Created cudaStream %ld for device #%d [%s]\n", StreamPool[s], DeviceList[i], DeviceProp.name);
#endif

	return 0;
}

// free-up all allocated resources
int CUDAEngine::Destroy()
{
	cudaDeviceSynchronize();

	// stop the profiler
//	if (cudaEngine->ProfilingMode == CUDAPROF_ALL) cudaEngine->StopProfiler();

	for (map<string, cudaEvent>::iterator i = EventPool.begin(); i != EventPool.end(); i++)
		cudaEventDestroy(i->second.event);

	for (int i = 0; i < StreamPool.size(); i++) cudaStreamDestroy(StreamPool[i]);

	map<void *, cudaPointerAttr >::iterator p;
	while (MemTable.begin() != MemTable.end()) Free((MemTable.begin())->first);

	StreamPool.clear();
	MemTable.clear();
	KernelTable.clear();
	EventPool.clear();
	TexturePtrs.clear();

	return 0;
}

// enable the CUDA hardware performance counter
void CUDAEngine::StartProfiler()
{
	cudaDeviceSynchronize();
	cudaProfilerStart();
	fprintf(stderr, "[CUDA] Profiler started\n");
}

// disable the CUDA hardware performance counter
void CUDAEngine::StopProfiler()
{
	cudaDeviceSynchronize();
	cudaProfilerStop();
	fprintf(stderr, "[CUDA] Profiler stopped\n");
}

// check if a given pointer matches the pointer binded to a given texture
void* CUDAEngine::TexturePtr(string tag, void *ptr)
{
	if (ptr == NULL) // query mode
	{
		map<string, void *>::iterator p = TexturePtrs.find(tag);
		if (p == TexturePtrs.end()) return NULL;
		else return p->second;
	}
	else
	{
		TexturePtrs[tag] = ptr;
		return ptr;
	}
}

// return a cuda event
cudaEvent_t CUDAEngine::Event(string tag, unsigned int flags)
{
	map<string, cudaEvent>::iterator p = EventPool.find(tag);
	if (p == EventPool.end())
	{
		cudaEvent_t e;
		cudaEventCreateWithFlags(&e, flags);
		EventPool[tag] = cudaEvent(e);
		p = EventPool.find(tag);
	}
	return p->second.event;
}

// free allocated memory
// thanks to the bookkeeping it can handle all sorts of memory
void CUDAEngine::Free(void *ptr)
{
	map<void *, cudaPointerAttr >::iterator p;
	p = MemTable.find(ptr);
	if (p == MemTable.end()) return;
	if (p->second.memoryType == cudaMemoryTypeDevice)
	{
		cudaFree(ptr);
#ifdef _RAINN_LOG_MEM
		fprintf(stderr, "[CUDA] GPU[%d] - : %s\n", p->second.device, size2text(p->second.size).c_str());
#endif
	}
	else
	{
		cudaFreeHost(ptr);
#ifdef _RAINN_LOG_MEM
		fprintf(stderr, "[CUDA] HOST - : %s\n", size2text(p->second.size).c_str());
#endif
	}
	MemTable.erase(p);
}

// return size of an allocated buffer
size_t CUDAEngine::QueryMemSize(void *ptr)
{
	map<void *, cudaPointerAttr >::iterator p;
	p = MemTable.find(ptr);
	if (p == MemTable.end()) return 0;
	else return p->second.size;
}

// return pitch of an pitched buffer
size_t CUDAEngine::QueryMemPitch(void *ptr)
{
	map<void *, cudaPointerAttr >::iterator p;
	p = MemTable.find(ptr);
	if (p == MemTable.end()) return 0;
	else return p->second.pitch;
}

// return height of an pitched buffer
size_t CUDAEngine::QueryMemHeight(void *ptr)
{
	map<void *, cudaPointerAttr >::iterator p;
	p = MemTable.find(ptr);
	if (p == MemTable.end()) return 0;
	else return p->second.height;
}

// return max block size of a kernel
size_t CUDAEngine::QueryKernelBlockSize(void *ptr)
{
	map<void *, cudaFuncAttributes >::iterator p;
	p = KernelTable.find(ptr);
	if (p == KernelTable.end())
	{
		cudaFuncAttributes attr;
		cudaFuncGetAttributes(&attr, ptr);
		KernelTable[ptr] = attr;
		if (attr.maxThreadsPerBlock == 0) fprintf(stderr, "[CUDA] Querying block size returned 0, something might have crashed somewhere else.\n");
		return attr.maxThreadsPerBlock;
	}
	else return p->second.maxThreadsPerBlock;
}

string CUDAEngine::size2text(size_t size)
{
	string Unit;
	float Size;
	if (size < 1024)
	{
		Unit = "B";
		Size = size;
	}
	else if (size < 1024 * 1024)
	{
		Unit = "KB";
		Size = size / 1024.0;
	}
	else
	{
		Unit = "MB";
		Size = size / 1024.0 / 1024.0;
	}
	char str[256];
	sprintf(str, "%.1f %s", Size, Unit.c_str());
	return string(str);
}

/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_CUDA_ENGINE_H
#define PS_CUDA_ENGINE_H

#include "pdps_cuda.h"
#include "pointers.h"
#include "cuda_occupancy.h"

namespace PDPS_NS {

	struct cudaPointerAttr : public cudaPointerAttributes
	{
		size_t size;
		size_t pitch;
		size_t height;
		string tag;
		cudaPointerAttr()
		{
			size = 0;
			pitch = 0;
			height = 0;
		}
	};

	struct cudaEvent
	{
		int         lock;
		cudaEvent_t event;
		cudaEvent() : lock(0), event(0) {}
		cudaEvent(cudaEvent_t eve) : lock(0), event(eve) {}
	};

	// CUDA management layer
	class CUDAEngine : public Pointers
	{
		friend class CUDAOccupancy;
	public:
		vector<cudaDeviceProp>           DeviceProperty;
		map<void *, cudaPointerAttr>     MemTable;
		map<void *, cudaFuncAttributes>  KernelTable;
		map<string, cudaEvent>           EventPool;
		map<string, void *>              TexturePtrs;
		string size2text(size_t size);

	public:
		int                              ProfilingMode;
		CUDAOccupancy                    OccuCalc;
		vector<int>				         DevicePool;
		vector<cudaStream_t>	         StreamPool;

		CUDAEngine(class PDPS *, int Device, string Profile);
		~CUDAEngine();
		int    Init(int Device, string Profile);
		int    Destroy();
		void   Free(void *ptr);
		void   StartProfiler();
		void   StopProfiler();
		size_t QueryMemSize(void *ptr);
		size_t QueryMemPitch(void *ptr);
		size_t QueryMemHeight(void *ptr);
		size_t QueryKernelBlockSize(void *ptr);
		void   Occupancy(int device_id, void *kernel, size_t dynamic_shmem, enum cudaFuncCache shmemConfig);
		cudaEvent_t Event(string tag, unsigned int flags = cudaEventDefault);
		void*  TexturePtr(string tag, void *ptr = NULL);

		// allocate 1d memory on device
		template <typename TYPE> TYPE *MallocDevice(string tag, size_t nElem)
		{
			TYPE *ptr;
			if (cudaMalloc(&ptr, nElem * sizeof(TYPE)) != cudaSuccess)
			{
				fprintf(stderr, "[CUDA] %s size = %s\n",
					cudaGetErrorString(cudaPeekAtLastError()),
					size2text(nElem * sizeof(TYPE)).c_str());
				cudaDeviceReset();
				exit(0);
			}
			cudaPointerAttr attr;
			cudaPointerGetAttributes(&attr, ptr);
			attr.size = nElem * sizeof(TYPE);
			attr.pitch = attr.size;
			attr.height = 1;
			attr.tag = tag;
			MemTable[ptr] = attr;
#ifdef _RAINN_LOG_MEM
			fprintf(stderr, "[CUDA] GPU[%d] + : %s\n", attr.device, size2text(attr.size).c_str());
#endif
			return ptr;
		}

		// allocate pitches quasi-2d memory on device
		template <typename TYPE> TYPE *MallocDevicePitch(string tag, size_t &pitch, size_t width, size_t height)
		{
			TYPE *ptr;
			if (cudaMallocPitch(&ptr, &pitch, width * sizeof(TYPE), height) != cudaSuccess)
			{
				fprintf(stderr, "[CUDA] %s\n", cudaGetErrorString(cudaPeekAtLastError()));
				cudaDeviceReset();
				exit(0);
			}
			cudaPointerAttr attr;
			cudaPointerGetAttributes(&attr, ptr);
			attr.pitch = pitch;
			attr.height = height;
			attr.size = pitch * height;
			attr.tag = tag;
			MemTable[ptr] = attr;
#ifdef _RAINN_LOG_MEM
			fprintf(stderr, "[CUDA] GPU[%d] + : %s\n", attr.device, size2text(attr.size).c_str());
#endif
			return ptr;
		}

		// resize allocated memory on device
		// can either grow/reduce
		// option: copy the original content
		// option: zero-out the new space
		template <typename TYPE> void ReallocDevice(string tag, TYPE*& ptr, size_t nElem, bool copy = true, bool zero = false)
		{
			map<void *, cudaPointerAttr >::iterator p;
			p = MemTable.find(ptr);
			if (p == MemTable.end()) // if NULL pointer
			{
				ptr = MallocDevice<TYPE>(tag, nElem);
				if (zero) cudaMemset(ptr, 0, nElem * sizeof(TYPE));
				return;
			}
			else if (p->second.memoryType == cudaMemoryTypeHost)
			{
				fprintf(stderr, "[CUDA] Host memory cannot be re-allocated to Device memory %p\n", ptr);
				cudaDeviceReset();
				exit(0);
			}
			else if (p->second.size == nElem * sizeof(TYPE)) return;

			if (copy)
			{
				TYPE* ptr_new = MallocDevice<TYPE>(tag, nElem);
				int   new_len = nElem * sizeof(TYPE);
				int   old_len = p->second.size;
				cudaMemcpy(ptr_new, ptr, min(old_len, new_len), cudaMemcpyDefault);
				Free(ptr);
				ptr = ptr_new;
			}
			else
			{
				Free(ptr);
				ptr = MallocDevice<TYPE>(tag, nElem);
				if (zero) cudaMemset(ptr, 0, nElem * sizeof(TYPE));
			}
		}

		// resize pitched memory on device
		template <typename TYPE> void ReallocDevicePitch(string tag, TYPE*& ptr, size_t &pitch, size_t width, size_t height, bool copy = true, bool zero = false)
		{
			map<void *, cudaPointerAttr >::iterator p;
			p = MemTable.find(ptr);
			if (p == MemTable.end()) // if NULL pointer
			{
				ptr = MallocDevicePitch<TYPE>(tag, pitch, width, height);
				if (zero) cudaMemset(ptr, 0, pitch * height);
				return;
			}
			else if (p->second.memoryType == cudaMemoryTypeHost)
			{
				fprintf(stderr, "[CUDA] Host memory cannot be re-allocated to Device memory %p\n", ptr);
				cudaDeviceReset();
				exit(0);
			}
			else if (p->second.pitch == width * sizeof(TYPE) && p->second.height == height) return;

			if (copy)
			{
				TYPE* ptr_new = MallocDevicePitch<TYPE>(tag, pitch, width, height);
				int   new_width = pitch;
				int   old_width = p->second.pitch;
				int   new_height = height;
				int   old_height = p->second.height;
				cudaMemcpy2D(ptr_new, pitch, ptr, p->second.pitch, min(old_width, new_width), min(old_height, new_height), cudaMemcpyDefault);
				Free(ptr);
				ptr = ptr_new;
			}
			else
			{
				Free(ptr);
				ptr = MallocDevicePitch<TYPE>(tag, pitch, width, height);
				if (zero) cudaMemset(ptr, 0, pitch * height * sizeof(TYPE));
			}
		}

		// allocate 1d host-pinned memory
		template <typename TYPE> TYPE *MallocHost(string tag, size_t nElem)
		{
			TYPE *ptr;
			if (cudaMallocHost(&ptr, nElem * sizeof(TYPE), cudaHostAllocMapped) != cudaSuccess)
			{
				fprintf(stderr, "[CUDA] %s size = %s\n",
					cudaGetErrorString(cudaPeekAtLastError()),
					size2text(nElem * sizeof(TYPE)).c_str());
				cudaDeviceReset();
				exit(0);
			}
			cudaPointerAttr attr;
			cudaPointerGetAttributes(&attr, ptr);
			attr.size = nElem * sizeof(TYPE);
			attr.pitch = attr.size;
			attr.height = 1;
			attr.tag = tag;
			MemTable[ptr] = attr;
#ifdef _RAINN_LOG_MEM
			fprintf(stderr, "[CUDA] HOST + : %s\n", size2text(attr.size).c_str());
#endif
			return ptr;
		}

		// allocate quasi-2d host pinned memory
		template <typename TYPE> TYPE *MallocHostPitch(string tag, size_t &pitch, size_t width, size_t height)
		{
			TYPE *ptr;
			pitch = ceiling(width * sizeof(TYPE)+127, 128);
			if (cudaMallocHost(&ptr, pitch * height, cudaHostAllocMapped) != cudaSuccess)
			{
				fprintf(stderr, "[CUDA] %s size = %s\n",
					cudaGetErrorString(cudaPeekAtLastError()),
					size2text(pitch * height).c_str());
				cudaDeviceReset();
				exit(0);
			}
			cudaPointerAttr attr;
			cudaPointerGetAttributes(&attr, ptr);
			attr.size = pitch * height;
			attr.pitch = pitch;
			attr.height = height;
			attr.tag = tag;
			MemTable[ptr] = attr;
#ifdef _RAINN_LOG_MEM
			fprintf(stderr, "[CUDA] HOST + : %s\n", size2text(attr.size).c_str());
#endif
			return ptr;
		}

		// reallocate host pinned memory
		template <typename TYPE> void ReallocHost(string tag, TYPE*& ptr, size_t nElem, bool copy = true, bool zero = false)
		{
			map<void *, cudaPointerAttr >::iterator p;
			p = MemTable.find(ptr);
			if (p == MemTable.end()) // if NULL pointer
			{
				ptr = MallocHost<TYPE>(tag, nElem);
				if (zero) memset(ptr, 0, nElem * sizeof(TYPE));
				return;
			}
			else if (p->second.memoryType == cudaMemoryTypeDevice)
			{
				fprintf(stderr, "[CUDA] Device memory cannot be re-allocated to Host memory %p\n", ptr);
				cudaDeviceReset();
				exit(0);
			}
			else if (p->second.size == nElem * sizeof(TYPE)) return;

			if (copy)
			{
				TYPE* ptr_new = MallocHost<TYPE>(tag, nElem);
				int   new_len = nElem * sizeof(TYPE);
				int   old_len = p->second.size;
				memcpy(ptr_new, ptr, min(old_len, new_len));
				Free(ptr);
				ptr = ptr_new;
			}
			else
			{
				Free(ptr);
				ptr = MallocHost<TYPE>(tag, nElem);
				if (zero) memset(ptr, 0, nElem * sizeof(TYPE));
			}
		}

		// pin existing host memory, allocate device buffers
		template<class TYPE>
		TYPE* PinHost(TYPE *&devHostArray, TYPE *ptrHostArray, size_t size)
		{
			if (devHostArray)
			{
				fprintf(stderr, "[CUDA] Previously pinned memory must be released before pinning new memory\n");
				return NULL;
			}
			else if (!ptrHostArray)
			{
				fprintf(stderr, "[CUDA] host memory submitted for pinning is null.\n");
				return NULL;
			}

			cudaHostRegister(ptrHostArray, size, cudaHostRegisterMapped);
			//		cudaHostGetDevicePointer( &devHostArray, ptrHostArray, 0 );
			devHostArray = MallocDevice<TYPE>("::HostArray", size / sizeof(TYPE));

#ifdef _RAINN_LOG_MEM
			fprintf(stderr, "[CUDA] Per-atom array pinned %p size = %d\n", ptrHostArray, size);
#endif

			return ptrHostArray;
		}

		// release pinned host memory
		template<class TYPE>
		void UnpinHost(TYPE *ptrHostArray, TYPE *&ptrHostArrayOld, TYPE *&devHostArray)
		{
			if (!ptrHostArray || !ptrHostArrayOld)
			{
				fprintf(stderr, "[CUDA] host memory submitted for unpinning is null.\n");
				return;
			}
			if (ptrHostArray && ptrHostArray != ptrHostArrayOld)
			{
				fprintf(stderr, "[CUDA] Previously pinned memory reallocated before been released\n");
				fprintf(stderr, "[CUDA] \tptrHostArray registered: %p, now: %p\n", ptrHostArrayOld, ptrHostArray);
				return;
			}

			cudaHostUnregister(ptrHostArray);
			Free(devHostArray);
			devHostArray = NULL;
			ptrHostArrayOld = NULL;

#ifdef _RAINN_LOG_MEM
			fprintf(stderr, "[CUDA] Per-atom array %p unpinned\n", ptrHostArray);
#endif
		}

		// return current time in nanosecond precision
		double GetTime()
		{
			//struct timespec time;
			//clock_gettime(CLOCK_REALTIME, &time);
			//return (double)time.tv_sec + (double)time.tv_nsec * 1.0e-9;
			clock_t current = clock();
			return (double)current;
		}
	};

}

#endif

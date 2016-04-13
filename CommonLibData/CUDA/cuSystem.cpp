/***********************************************************************
 *  File:       cuSystem.cpp
 *
 *  Purpose:    Implementation of a a CUDA-info class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuSystem.hpp"

namespace gem {

int ConvertSMVer2Cores(int major, int minor)
{
    // Defines for GPU Architecture types (using the SM version to
    // determine the # of cores per SM
    typedef struct {
        int SM;     // 0xMm (hexidecimal notation),
                    // M = SM Major version,
                    // and m = SM minor version
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] =
    { { 0x10,  8 },
      { 0x11,  8 },
      { 0x12,  8 },
      { 0x13,  8 },
      { 0x20, 32 },
      { 0x21, 48 },
      {   -1, -1 }
    };

    int index = 0;
    while (nGpuArchCoresPerSM[index].SM != -1) {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor) ) {
            return nGpuArchCoresPerSM[index].Cores;
        }
        index++;
    }
    printf("MapSMtoCores undefined SMversion %d.%d!\n", major, minor);
    return -1;
}

cuSystem::cuSystem()
{
    CUDA_SAFE_CALL(cudaGetDeviceCount(&_numDevice));

    if (_numDevice > 0) {
        _deviceProp = new cudaDeviceProp [_numDevice];

        for (int i = 0; i < _numDevice; i++) {
            CUDA_SAFE_CALL(cudaGetDeviceProperties(&_deviceProp[i],i));
        }

        CUDA_SAFE_CALL(cudaDriverGetVersion (&_verDriver ));
        CUDA_SAFE_CALL(cudaRuntimeGetVersion(&_verRunTime));
    }
    else {
        _deviceProp = NULL;
    }
}

cuSystem::~cuSystem()
{
    if (_deviceProp != NULL) {
        delete [] _deviceProp;
    }
}

void cuSystem::resetDevice(void)
{
    cudaDeviceReset();
}

void cuSystem::printDevInfo()
{
    if (_numDevice == 0) {
        printf("There is no device supporting CUDA \n");
        return;
    }
    else {
        printf("Found %d CUDA Capable device(s)\n", _numDevice);
        printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
            _verDriver /1000, (_verDriver %100)/10,
            _verRunTime/1000, (_verRunTime%100)/10);
    }

    for (int i = 0; i < _numDevice; i++) {
        printf("\nDevice %d: \"%s\"\n", i, _deviceProp[i].name);
        printf("  CUDA Capability Major/Minor version number:    %d.%d\n",
            _deviceProp[i].major,
            _deviceProp[i].minor);
        printf("  Total amount of global memory:                 %.0f MBytes (%llu bytes)\n",
            (float) _deviceProp[i].totalGlobalMem/1048576.0f,
            (unsigned long long) _deviceProp[i].totalGlobalMem);
        printf("  (%2d) Multiprocessors x (%2d) CUDA Cores/MP:     %d CUDA Cores\n",
            _deviceProp[i].multiProcessorCount,
            ConvertSMVer2Cores(_deviceProp[i].major, _deviceProp[i].minor),
            ConvertSMVer2Cores(_deviceProp[i].major, _deviceProp[i].minor) *
            _deviceProp[i].multiProcessorCount);
        printf("  GPU Clock Speed:                               %.2f GHz\n",
            (float) _deviceProp[i].clockRate * 1e-6f);
        printf("  Memory Clock rate:                             %.2f Mhz\n",
            (float) _deviceProp[i].memoryClockRate * 1e-3f);
        printf("  Memory Bus Width:                              %d-bit\n",
            _deviceProp[i].memoryBusWidth);
        printf("  L2 Cache Size:                                 %d bytes\n",
            _deviceProp[i].l2CacheSize);
        printf("  Max Texture Dimension Size (x,y,z)             1D=(%d), 2D=(%d,%d), 3D=(%d,%d,%d)\n",
            _deviceProp[i].maxTexture1D,
            _deviceProp[i].maxTexture2D[0],
            _deviceProp[i].maxTexture2D[1],
            _deviceProp[i].maxTexture3D[0],
            _deviceProp[i].maxTexture3D[1],
            _deviceProp[i].maxTexture3D[2]);
        printf("  Max Layered Texture Size (dim) x layers        1D=(%d) x %d, 2D=(%d,%d) x %d\n",
            _deviceProp[i].maxTexture1DLayered[0], _deviceProp[i].maxTexture1DLayered[1],
            _deviceProp[i].maxTexture2DLayered[0], _deviceProp[i].maxTexture2DLayered[1],
            _deviceProp[i].maxTexture2DLayered[2]);
        printf("  Total amount of constant memory:               %llu bytes\n",
            (unsigned long long) _deviceProp[i].totalConstMem);
        printf("  Total amount of shared memory per block:       %llu bytes\n",
            (unsigned long long) _deviceProp[i].sharedMemPerBlock);
        printf("  Total number of registers available per block: %d\n",
            _deviceProp[i].regsPerBlock);
        printf("  Warp size:                                     %d\n",
            _deviceProp[i].warpSize);
        printf("  Maximum number of threads per block:           %d\n",
            _deviceProp[i].maxThreadsPerBlock);
        printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
            _deviceProp[i].maxThreadsDim[0],
            _deviceProp[i].maxThreadsDim[1],
            _deviceProp[i].maxThreadsDim[2]);
        printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
            _deviceProp[i].maxGridSize[0],
            _deviceProp[i].maxGridSize[1],
            _deviceProp[i].maxGridSize[2]);
        printf("  Maximum memory pitch:                          %llu bytes\n",
            (unsigned long long) _deviceProp[i].memPitch);
        printf("  Texture alignment:                             %llu bytes\n",
            (unsigned long long) _deviceProp[i].textureAlignment);
        printf("  Concurrent copy and execution:                 %s with %d copy engine(s)\n",
            (_deviceProp[i].deviceOverlap ? "Yes" : "No"),
            _deviceProp[i].asyncEngineCount);
        printf("  Run time limit on kernels:                     %s\n",
            _deviceProp[i].kernelExecTimeoutEnabled ? "Yes" : "No");
        printf("  Integrated GPU sharing Host Memory:            %s\n",
            _deviceProp[i].integrated ? "Yes" : "No");
        printf("  Support host page-locked memory mapping:       %s\n",
            _deviceProp[i].canMapHostMemory ? "Yes" : "No");
        printf("  Concurrent kernel execution:                   %s\n",
            _deviceProp[i].concurrentKernels ? "Yes" : "No");
        printf("  Alignment requirement for Surfaces:            %s\n",
            _deviceProp[i].surfaceAlignment ? "Yes" : "No");
        printf("  Device has ECC support enabled:                %s\n",
            _deviceProp[i].ECCEnabled ? "Yes" : "No");
        printf("  Device is using TCC driver mode:               %s\n",
            _deviceProp[i].tccDriver ? "Yes" : "No");
        printf("  Device supports Unified Addressing (UVA):      %s\n",
            _deviceProp[i].unifiedAddressing ? "Yes" : "No");
        printf("  Device PCI Bus ID / PCI location ID:           %d / %d\n",
            _deviceProp[i].pciBusID, _deviceProp[i].pciDeviceID );
        const char *sComputeMode[] = {
            "Default (multiple host threads can use ::cudaSetDevice() with device simultaneously)",
            "Exclusive (only one host thread in one process is able to use ::cudaSetDevice() with this device)",
            "Prohibited (no host thread can use ::cudaSetDevice() with this device)",
            "Exclusive Process (many threads in one process is able to use ::cudaSetDevice() with this device)",
            "Unknown",
            NULL
        };
        printf("  Compute Mode:\n");
        printf("     < %s >\n\n", sComputeMode[_deviceProp[i].computeMode]);
    }
}

void cuSystem::printMemInfo(int devID)
{
    if (devID < _numDevice) {
        size_t memFree, memTotal;

        CUDA_SAFE_CALL(cudaSetDevice(devID));
        CUDA_SAFE_CALL(cudaMemGetInfo(&memFree, &memTotal));

        std::cout << "Device " << devID << " GPU memory usage: \n"
                  << "   total = " << (double) memTotal/1048576.0
                  << " MB" << std::endl
                  << "    free = " << (double) memFree/1048576.0
                  << " MB" << std::endl << std::endl;
    }
    else {
        std::cout << "The input devID = " << devID << " is too large!\n"
                  << std::endl;
    }
}

void cuSystem::printThroughput(size_t numElements, float time)
{
    float   throughput = 2.0f * sizeof(float) * (float) numElements / (1024*1024*1024) / time;
    printf("Throughput = %.4f GB/s\n\n", throughput);
}

} // namespace gem

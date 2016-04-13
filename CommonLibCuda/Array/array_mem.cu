/***********************************************************************
 *  File:       array_mem.cu
 *
 *  Purpose:    Implementation of array-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "array.cuh"

#include <curand_kernel.h>

#ifndef NDEBUG

/*#ifdef __GEM_USE_STD11__
    #include <mutex>

namespace gem {

    std::mutex      mutexAllocCUDA;

} // namespace gem

#else*/
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wconversion"
    #include "boost/thread.hpp"
    #pragma GCC diagnostic pop

namespace gem {

    boost::mutex    mutexAllocCUDA;

} // namespace gem

#endif

namespace gem {

static int      countAllocDev  = 0;
static int      countAllocArr  = 0;
static int      countAllocHost = 0;

} // namespace gem

//#endif

namespace gem {

/*****************************************
 * memcpy / memset / memcheck
 ****************************************/

// memcpy device -> device
template <typename T>
void cuda_array_memcpy_d2d(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    CUDA_SAFE_CALL( cudaMemcpy(arrayDst, arraySrc, nRow*sizeof(T),
                               cudaMemcpyDeviceToDevice) );
}

// instantiation
template
void cuda_array_memcpy_d2d<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_memcpy_d2d<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_memcpy_d2d<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_memcpy_d2d<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_memcpy_d2d<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_memcpy_d2d<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// memcpy device -> host
template <typename T>
void cuda_array_memcpy_d2h(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    CUDA_SAFE_CALL( cudaMemcpy(arrayDst, arraySrc, nRow*sizeof(T),
                               cudaMemcpyDeviceToHost) );
}

// instantiation
template
void cuda_array_memcpy_d2h<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_memcpy_d2h<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_memcpy_d2h<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_memcpy_d2h<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_memcpy_d2h<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_memcpy_d2h<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// memcpy host -> device
template <typename T>
void cuda_array_memcpy_h2d(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    CUDA_SAFE_CALL( cudaMemcpy(arrayDst, arraySrc, nRow*sizeof(T),
                               cudaMemcpyHostToDevice) );
}

// instantiation
template
void cuda_array_memcpy_h2d<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_memcpy_h2d<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_memcpy_h2d<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_memcpy_h2d<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_memcpy_h2d<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_memcpy_h2d<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// memcpy host -> host
template <typename T>
void cuda_array_memcpy_h2h(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    CUDA_SAFE_CALL( cudaMemcpy(arrayDst, arraySrc, nRow*sizeof(T),
                               cudaMemcpyHostToHost) );
}

// instantiation
template
void cuda_array_memcpy_h2h<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_memcpy_h2h<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_memcpy_h2h<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_memcpy_h2h<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_memcpy_h2h<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_memcpy_h2h<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// memcpy 2D device -> device's cudaArray
void cuda_array_memcpy_d2d(cudaArray* const arrayDst, const float* const arraySrc, size_t nRow, size_t nCol)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);

    CUDA_SAFE_CALL( cudaMemcpyToArray(arrayDst, 0, 0, arraySrc,
                        nRow*nCol*sizeof(float),
                        cudaMemcpyDeviceToDevice) );
}

// memcpy 2D host -> device's cudaArray
void cuda_array_memcpy_h2d(cudaArray* const arrayDst, const float* const arraySrc, size_t nRow, size_t nCol)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);

    CUDA_SAFE_CALL( cudaMemcpyToArray(arrayDst, 0, 0, arraySrc,
                        nRow*nCol*sizeof(float),
                        cudaMemcpyHostToDevice) );
}

// memcpy 3D device -> device's cudaArray
void cuda_array_memcpy_d2d(cudaArray* const arrayDst, const float* const arraySrc, size_t nRow, size_t nCol, size_t nSec)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    cudaExtent           volumeSize = make_cudaExtent(nSec, nCol, nRow);
    cudaMemcpy3DParms    copyParams;

    copyParams.srcArray = 0;
    copyParams.srcPos   = make_cudaPos(0, 0, 0);
    copyParams.srcPtr   = make_cudaPitchedPtr((void *) (intptr_t) arraySrc,
                                              volumeSize.width*sizeof(float),
                                              volumeSize.width,
                                              volumeSize.height);
    copyParams.dstArray = arrayDst;
    copyParams.dstPos   = make_cudaPos(0, 0, 0);
    copyParams.dstPtr   = make_cudaPitchedPtr(0, 0, 0, 0);
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyDeviceToDevice;

    CUDA_SAFE_CALL( cudaMemcpy3D(&copyParams) );
}

// memcpy 3D host -> device's cudaArray
void cuda_array_memcpy_h2d(cudaArray* const arrayDst, const float* const arraySrc, size_t nRow, size_t nCol, size_t nSec)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    cudaExtent           volumeSize = make_cudaExtent(nSec, nCol, nRow);
    cudaMemcpy3DParms    copyParams;

    copyParams.srcArray = 0;
    copyParams.srcPos   = make_cudaPos(0, 0, 0);
    copyParams.srcPtr   = make_cudaPitchedPtr((void *) (intptr_t) arraySrc,
                                              volumeSize.width*sizeof(float),
                                              volumeSize.width,
                                              volumeSize.height);
    copyParams.dstArray = arrayDst;
    copyParams.dstPos   = make_cudaPos(0, 0, 0);
    copyParams.dstPtr   = make_cudaPitchedPtr(0, 0, 0, 0);
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;

    CUDA_SAFE_CALL( cudaMemcpy3D(&copyParams) );
}

// memset
template <typename T>
void cuda_array_memset(T* const array, int value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    CUDA_SAFE_CALL( cudaMemset(array, value, nRow*sizeof(T)) );
}

// instantiation
template
void cuda_array_memset<int32_t        >(int32_t*         const array, int value, size_t nRow);
template
void cuda_array_memset<uint32_t       >(uint32_t*        const array, int value, size_t nRow);
template
void cuda_array_memset<int64_t        >(int64_t*         const array, int value, size_t nRow);
template
void cuda_array_memset<uint64_t       >(uint64_t*        const array, int value, size_t nRow);
template
void cuda_array_memset<float          >(float*           const array, int value, size_t nRow);
template
void cuda_array_memset<double         >(double*          const array, int value, size_t nRow);
template
void cuda_array_memset<cuFloatComplex >(cuFloatComplex*  const array, int value, size_t nRow);
template
void cuda_array_memset<cuDoubleComplex>(cuDoubleComplex* const array, int value, size_t nRow);

// check memory leak
void cuda_array_memcheck(void)
{
    require(countAllocDev  == 0, "array_memcheck: countAllocDev = " + num2str(countAllocDev));
    require(countAllocArr  == 0, "array_memcheck: countAllocArr = " + num2str(countAllocArr));
    require(countAllocHost == 0, "array_memcheck: countAllocHost = " + num2str(countAllocHost));
}

/*****************************************
 * Array allocation/deallocation
 *       w/o zero initialization
 ****************************************/

// device memory: cuArray
void cuda_arrayArr_new(cudaArray*& array, size_t nRow, size_t nCol, size_t nSec)
{
    assert(array == NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    cudaChannelFormatDesc    channelDesc = cudaCreateChannelDesc<float>();
    cudaExtent               volumeSize  = make_cudaExtent(nSec, nCol, nRow);

    CUDA_SAFE_CALL(cudaMalloc3DArray(&array, &channelDesc, volumeSize));

#ifndef NDEBUG
mutexAllocCUDA.lock();
        countAllocArr++;
mutexAllocCUDA.unlock();
#endif
}

void cuda_arrayArr_new(cudaArray*& array, size_t nRow, size_t nCol)
{
    assert(array == NULL);
    assert(nRow > 0 && nCol > 0);

    cudaChannelFormatDesc    channelDesc = cudaCreateChannelDesc<float>();

    CUDA_SAFE_CALL( cudaMallocArray(&array, &channelDesc, nCol, nRow) );

#ifndef NDEBUG
mutexAllocCUDA.lock();
        countAllocArr++;
mutexAllocCUDA.unlock();
#endif
}

void cuda_arrayArr_delete(cudaArray*& array)
{
    if (array != NULL) {
        CUDA_SAFE_CALL( cudaFreeArray(array) );
        array = NULL;

#ifndef NDEBUG
mutexAllocCUDA.lock();
        countAllocArr--;
mutexAllocCUDA.unlock();
#endif
    }
    else {
#ifndef NDEBUG
        WARNING("cuda_arrayArr_delete", "attempt to delete a NULL object");
#endif
    }
}

// size(array1D) = [nRow]
template <typename T>
void cuda_arrayDev_new(T*& array1D, size_t nRow)
{
    assert(array1D == NULL);
    assert(nRow > 0);

    CUDA_SAFE_CALL( cudaMalloc((void **) &array1D, nRow*sizeof(T)) );

#ifndef NDEBUG
mutexAllocCUDA.lock();
    countAllocDev++;
mutexAllocCUDA.unlock();
#endif
}

// instantiation
template
void cuda_arrayDev_new<int32_t        >(int32_t*&         array1D, size_t nRow);
template
void cuda_arrayDev_new<uint32_t       >(uint32_t*&        array1D, size_t nRow);
template
void cuda_arrayDev_new<int64_t        >(int64_t*&         array1D, size_t nRow);
template
void cuda_arrayDev_new<uint64_t       >(uint64_t*&        array1D, size_t nRow);
template
void cuda_arrayDev_new<float          >(float*&           array1D, size_t nRow);
template
void cuda_arrayDev_new<double         >(double*&          array1D, size_t nRow);
template
void cuda_arrayDev_new<cuFloatComplex >(cuFloatComplex*&  array1D, size_t nRow);
template
void cuda_arrayDev_new<cuDoubleComplex>(cuDoubleComplex*& array1D, size_t nRow);
template
void cuda_arrayDev_new<curandState    >(curandState*&     array1D, size_t nRow);
template
void cuda_arrayDev_new<float*         >(float**&          array1D, size_t nRow);
template
void cuda_arrayDev_new<double*        >(double**&         array1D, size_t nRow);

// size(array1D) = [nRow]
template <typename T>
void cuda_arrayDev_new_zero(T*& array1D, size_t nRow)
{
    cuda_arrayDev_new(array1D, nRow);
    cuda_array_memset(array1D, 0, nRow);
}

// instantiation
template
void cuda_arrayDev_new_zero<int32_t        >(int32_t*&         array1D, size_t nRow);
template
void cuda_arrayDev_new_zero<uint32_t       >(uint32_t*&        array1D, size_t nRow);
template
void cuda_arrayDev_new_zero<int64_t        >(int64_t*&         array1D, size_t nRow);
template
void cuda_arrayDev_new_zero<uint64_t       >(uint64_t*&        array1D, size_t nRow);
template
void cuda_arrayDev_new_zero<float          >(float*&           array1D, size_t nRow);
template
void cuda_arrayDev_new_zero<double         >(double*&          array1D, size_t nRow);
template
void cuda_arrayDev_new_zero<cuFloatComplex >(cuFloatComplex*&  array1D, size_t nRow);
template
void cuda_arrayDev_new_zero<cuDoubleComplex>(cuDoubleComplex*& array1D, size_t nRow);
template
void cuda_arrayDev_new_zero<curandState    >(curandState*&     array1D, size_t nRow);
template
void cuda_arrayDev_new_zero<float*         >(float**&          array1D, size_t nRow);
template
void cuda_arrayDev_new_zero<double*        >(double**&         array1D, size_t nRow);

// size(array1D) = [nRow]
template <typename T>
void cuda_arrayDev_delete(T*& array1D)
{
    if (array1D != NULL) {
        CUDA_SAFE_CALL( cudaFree(array1D) );
        array1D = NULL;

#ifndef NDEBUG
mutexAllocCUDA.lock();
        countAllocDev--;
mutexAllocCUDA.unlock();
#endif
    }
    else {
#ifndef NDEBUG
        WARNING("cuda_arrayDev_delete", "attempt to delete a NULL object");
#endif
    }
}

// instantiation
template
void cuda_arrayDev_delete<int32_t        >(int32_t*&         array1D);
template
void cuda_arrayDev_delete<uint32_t       >(uint32_t*&        array1D);
template
void cuda_arrayDev_delete<int64_t        >(int64_t*&         array1D);
template
void cuda_arrayDev_delete<uint64_t       >(uint64_t*&        array1D);
template
void cuda_arrayDev_delete<float          >(float*&           array1D);
template
void cuda_arrayDev_delete<double         >(double*&          array1D);
template
void cuda_arrayDev_delete<cuFloatComplex >(cuFloatComplex*&  array1D);
template
void cuda_arrayDev_delete<cuDoubleComplex>(cuDoubleComplex*& array1D);
template
void cuda_arrayDev_delete<curandState    >(curandState*&     array1D);
template
void cuda_arrayDev_delete<float*         >(float**&          array1D);
template
void cuda_arrayDev_delete<double*        >(double**&         array1D);

// size(array2D) = [nRow nCol]
template <typename T>
void cuda_arrayDev_new(T**& array2D, size_t nRow, size_t nCol)
{
    assert(array2D == NULL);
    assert(nRow > 0);
    assert(nCol > 0);

    CUDA_SAFE_CALL( cudaMalloc((void **) &array2D, nRow*sizeof(*array2D)) );

    for (size_t i = 0; i < nRow; i++) {
        cuda_arrayDev_new(array2D[i], nCol);
    }

#ifndef NDEBUG
mutexAllocCUDA.lock();
    countAllocDev++;
mutexAllocCUDA.unlock();
#endif
}

// instantiation
template
void cuda_arrayDev_new<int32_t >(int32_t**&  array2D, size_t nRow, size_t nCol);
template
void cuda_arrayDev_new<uint32_t>(uint32_t**& array2D, size_t nRow, size_t nCol);
template
void cuda_arrayDev_new<int64_t >(int64_t**&  array2D, size_t nRow, size_t nCol);
template
void cuda_arrayDev_new<uint64_t>(uint64_t**& array2D, size_t nRow, size_t nCol);
template
void cuda_arrayDev_new<float   >(float**&    array2D, size_t nRow, size_t nCol);
template
void cuda_arrayDev_new<double  >(double**&   array2D, size_t nRow, size_t nCol);

// size(array2D) = [nRow nCol]
template <typename T>
void cuda_arrayDev_new_zero(T**& array2D, size_t nRow, size_t nCol)
{
    cuda_arrayDev_new_zero(array2D, nRow, nCol);

    for (size_t i = 0; i < nRow; i++) {
        cuda_array_memset(array2D[i], 0, nCol);
    }
}

// instantiation
template
void cuda_arrayDev_new_zero<int32_t >(int32_t**&  array2D, size_t nRow, size_t nCol);
template
void cuda_arrayDev_new_zero<uint32_t>(uint32_t**& array2D, size_t nRow, size_t nCol);
template
void cuda_arrayDev_new_zero<int64_t >(int64_t**&  array2D, size_t nRow, size_t nCol);
template
void cuda_arrayDev_new_zero<uint64_t>(uint64_t**& array2D, size_t nRow, size_t nCol);
template
void cuda_arrayDev_new_zero<float   >(float**&    array2D, size_t nRow, size_t nCol);
template
void cuda_arrayDev_new_zero<double  >(double**&   array2D, size_t nRow, size_t nCol);

// size(array2D) = [nRow nCol]
template <typename T>
void cuda_arrayDev_delete(T**& array2D, size_t nRow)
{
    assert(nRow > 0);

    if (array2D != NULL) {
        for (size_t i = 0; i < nRow; i++) {
            if (array2D[i] != NULL) {
                cuda_arrayDev_delete(array2D[i]);
            }
        }

        CUDA_SAFE_CALL( cudaFree(array2D) );
        array2D = NULL;

#ifndef NDEBUG
mutexAllocCUDA.lock();
        countAllocDev--;
mutexAllocCUDA.unlock();
#endif
    }
    else {
#ifndef NDEBUG
        WARNING("cuda_arrayDev_delete", "attempt to delete a NULL object");
#endif
    }
}

// instantiation
template
void cuda_arrayDev_delete<int32_t >(int32_t**&  array2D, size_t nRow);
template
void cuda_arrayDev_delete<uint32_t>(uint32_t**& array2D, size_t nRow);
template
void cuda_arrayDev_delete<int64_t >(int64_t**&  array2D, size_t nRow);
template
void cuda_arrayDev_delete<uint64_t>(uint64_t**& array2D, size_t nRow);
template
void cuda_arrayDev_delete<float   >(float**&    array2D, size_t nRow);
template
void cuda_arrayDev_delete<double  >(double**&   array2D, size_t nRow);

// size(array1D) = [nRow]
template <typename T>
void cuda_arrayHost_new(T*& array1D, size_t nRow)
{
    assert(array1D == NULL);
    assert(nRow > 0);

    CUDA_SAFE_CALL( cudaMallocHost((void **) &array1D, nRow*sizeof(T)) );

#ifndef NDEBUG
mutexAllocCUDA.lock();
    countAllocHost++;
mutexAllocCUDA.lock();
#endif
}

// instantiation
template
void cuda_arrayHost_new<int32_t >(int32_t*&  array1D, size_t nRow);
template
void cuda_arrayHost_new<uint32_t>(uint32_t*& array1D, size_t nRow);
template
void cuda_arrayHost_new<int64_t >(int64_t*&  array1D, size_t nRow);
template
void cuda_arrayHost_new<uint64_t>(uint64_t*& array1D, size_t nRow);
template
void cuda_arrayHost_new<float   >(float*&    array1D, size_t nRow);
template
void cuda_arrayHost_new<double  >(double*&   array1D, size_t nRow);

// size(array1D) = [nRow]
template <typename T>
void cuda_arrayHost_new_zero(T*& array1D, size_t nRow)
{
    cuda_arrayHost_new(array1D, nRow);
    array_memset(array1D, 0, nRow);
}

// instantiation
template
void cuda_arrayHost_new_zero<int32_t >(int32_t*&  array1D, size_t nRow);
template
void cuda_arrayHost_new_zero<uint32_t>(uint32_t*& array1D, size_t nRow);
template
void cuda_arrayHost_new_zero<int64_t >(int64_t*&  array1D, size_t nRow);
template
void cuda_arrayHost_new_zero<uint64_t>(uint64_t*& array1D, size_t nRow);
template
void cuda_arrayHost_new_zero<float   >(float*&    array1D, size_t nRow);
template
void cuda_arrayHost_new_zero<double  >(double*&   array1D, size_t nRow);

// size(array1D) = [nRow]
template <typename T>
void cuda_arrayHost_delete(T*& array1D)
{
    if (array1D != NULL) {
        CUDA_SAFE_CALL( cudaFreeHost(array1D) );
        array1D = NULL;

#ifndef NDEBUG
mutexAllocCUDA.lock();
        countAllocHost--;
mutexAllocCUDA.unlock();
#endif
    }
    else {
#ifndef NDEBUG
        WARNING("cuda_arrayHost_delete", "attempt to delete a NULL object");
#endif
    }
}

// instantiation
template
void cuda_arrayHost_delete<int32_t >(int32_t*&  array1D);
template
void cuda_arrayHost_delete<uint32_t>(uint32_t*& array1D);
template
void cuda_arrayHost_delete<int64_t >(int64_t*&  array1D);
template
void cuda_arrayHost_delete<uint64_t>(uint64_t*& array1D);
template
void cuda_arrayHost_delete<float   >(float*&    array1D);
template
void cuda_arrayHost_delete<double  >(double*&   array1D);

} // namespace gem

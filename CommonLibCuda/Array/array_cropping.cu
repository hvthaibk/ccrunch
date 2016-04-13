/***********************************************************************
 *  File:       array_cropping.cu
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

namespace gem {

/*****************************************
 * Cropping
 ****************************************/

// 1D
template <typename T> __global__
void dev_array_crop(const T* const arraySrc, size_t nRowSrc,
                          T* const arrayDst, size_t nRowDst,
                                             size_t nRowOff)
{
    size_t    iRowDst = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowDst < nRowDst) {
        arrayDst[iRowDst] = arraySrc[iRowDst+nRowOff];
    }
}

template <typename T>
void cuda_array_crop(const T* const arraySrc, size_t nRowSrc,
                           T* const arrayDst, size_t nRowDst,
                                              size_t nRowOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0);
    assert(nRowDst > 0);
    assert(nRowDst + nRowOff <= nRowSrc);

    dev_array_crop<<<iDivUp(nRowDst, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arraySrc, nRowSrc,
         arrayDst, nRowDst,
                    nRowOff);
}

// instantiation
template
void cuda_array_crop<int32_t >(const int32_t*  const arraySrc, size_t nRowSrc,
                                     int32_t*  const arrayDst, size_t nRowDst,
                                                               size_t nRowOff);
template
void cuda_array_crop<uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc,
                                     uint32_t* const arrayDst, size_t nRowDst,
                                                               size_t nRowOff);
template
void cuda_array_crop<int64_t >(const int64_t*  const arraySrc, size_t nRowSrc,
                                     int64_t*  const arrayDst, size_t nRowDst,
                                                               size_t nRowOff);
template
void cuda_array_crop<uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc,
                                     uint64_t* const arrayDst, size_t nRowDst,
                                                               size_t nRowOff);
template
void cuda_array_crop<float   >(const float*    const arraySrc, size_t nRowSrc,
                                     float*    const arrayDst, size_t nRowDst,
                                                               size_t nRowOff);
template
void cuda_array_crop<double  >(const double*   const arraySrc, size_t nRowSrc,
                                     double*   const arrayDst, size_t nRowDst,
                                                               size_t nRowOff);

// 2D
template <typename T> __global__
void dev_array_crop(const T* const arraySrc, size_t nRowSrc, size_t nColSrc,
                          T* const arrayDst, size_t nRowDst, size_t nColDst,
                                             size_t nRowOff, size_t nColOff)
{
    size_t    iRowDst = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iColDst = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowDst < nRowDst && iColDst < nColDst) {
        arrayDst[iRowDst*nColDst+iColDst] = arraySrc[(iRowDst+nRowOff)*nColSrc+iColDst+nColOff];
    }
}

template <typename T>
void cuda_array_crop(const T* const arraySrc, size_t nRowSrc, size_t nColSrc,
                           T* const arrayDst, size_t nRowDst, size_t nColDst,
                                              size_t nRowOff, size_t nColOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0 && nColSrc > 0);
    assert(nRowDst > 0 && nColDst > 0);
    assert(nRowDst + nRowOff <= nRowSrc);
    assert(nColDst + nColOff <= nColSrc);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nColDst, dimBlock.x),
                     iDivUp(nRowDst, dimBlock.y));

    dev_array_crop<<<dimGrid, dimBlock>>>
        (arraySrc, nRowSrc, nColSrc,
         arrayDst, nRowDst, nColDst,
                   nRowOff, nColOff);
}

// instantiation
template
void cuda_array_crop<int32_t >(const int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                     int32_t*  const arrayDst, size_t nRowDst, size_t nColDst,
                                                               size_t nRowOff, size_t nColOff);
template
void cuda_array_crop<uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                     uint32_t* const arrayDst, size_t nRowDst, size_t nColDst,
                                                               size_t nRowOff, size_t nColOff);
template
void cuda_array_crop<int64_t >(const int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                     int64_t*  const arrayDst, size_t nRowDst, size_t nColDst,
                                                               size_t nRowOff, size_t nColOff);
template
void cuda_array_crop<uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                     uint64_t* const arrayDst, size_t nRowDst, size_t nColDst,
                                                               size_t nRowOff, size_t nColOff);
template
void cuda_array_crop<float   >(const float*    const arraySrc, size_t nRowSrc, size_t nColSrc,
                                     float*    const arrayDst, size_t nRowDst, size_t nColDst,
                                                               size_t nRowOff, size_t nColOff);
template
void cuda_array_crop<double  >(const double*   const arraySrc, size_t nRowSrc, size_t nColSrc,
                                     double*   const arrayDst, size_t nRowDst, size_t nColDst,
                                                               size_t nRowOff, size_t nColOff);

// 3D
template <typename T> __global__
void dev_array_crop(const T* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                          T* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                             size_t nRowOff, size_t nColOff, size_t nSecOff)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    size_t    iRowDst;
    size_t    iColDst = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSecDst = blockDim.x * blockIdx.x + threadIdx.x;

    if (iColDst < nColDst && iSecDst < nSecDst) {
        for (iRowDst = 0; iRowDst < nRowDst; iRowDst++) {
            arrayDst[(iRowDst*nColDst+iColDst)*nSecDst+iSecDst] =
            arraySrc[((iRowDst+nRowOff)*nColSrc+(iColDst+nColOff))*nSecSrc+iSecDst+nSecOff];
        }
    }
#else
    size_t    iRowDst = blockDim.z * blockIdx.z + threadIdx.z;
    size_t    iColDst = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSecDst = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowDst < nRowDst && iColDst < nColDst && iSecDst < nSecDst) {
        arrayDst[(iRowDst*nColDst+iColDst)*nSecDst+iSecDst] =
        arraySrc[((iRowDst+nRowOff)*nColSrc+(iColDst+nColOff))*nSecSrc+iSecDst+nSecOff];
    }
#endif
}

template <typename T>
void cuda_array_crop(const T* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                           T* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                              size_t nRowOff, size_t nColOff, size_t nSecOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);
    assert(nRowDst > 0 && nColDst > 0 && nSecDst > 0);
    assert(nRowDst + nRowOff <= nRowSrc);
    assert(nColDst + nColOff <= nColSrc);
    assert(nSecDst + nSecOff <= nSecSrc);

#ifdef __GEM_CUDA_ARCH_HOST_130__
    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nSecDst, dimBlock.x),
                     iDivUp(nColDst, dimBlock.y));
#else
    dim3    dimBlock(BLOCK_3D_NSEC, BLOCK_3D_NCOL, BLOCK_3D_NROW);
    dim3    dimGrid (iDivUp(nSecDst, dimBlock.x),
                     iDivUp(nColDst, dimBlock.y),
                     iDivUp(nRowDst, dimBlock.z));
#endif

    dev_array_crop<<<dimGrid, dimBlock>>>
        (arraySrc, nRowSrc, nColSrc, nSecSrc,
         arrayDst, nRowDst, nColDst, nSecDst,
                   nRowOff, nColOff, nSecOff);
}

// instantiation
template
void cuda_array_crop<int32_t >(const int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                     int32_t*  const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                               size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_crop<uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                     uint32_t* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                               size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_crop<int64_t >(const int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                     int64_t*  const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                               size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_crop<uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                     uint64_t* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                               size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_crop<float   >(const float*    const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                     float*    const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                               size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_crop<double  >(const double*   const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                     double*   const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                               size_t nRowOff, size_t nColOff, size_t nSecOff);

} // namespace gem

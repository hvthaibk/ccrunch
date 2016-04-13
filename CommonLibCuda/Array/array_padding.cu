/***********************************************************************
 *  File:       array_padding.cu
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
 * Padding
 ****************************************/

// 1D
template <typename T> __global__
void dev_array_pad(const T* const arraySrc, size_t nRowSrc,
                         T* const arrayDst, size_t nRowDst,
                                            size_t nRowOff)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRowSrc) {
        arrayDst[iRow+nRowOff] = arraySrc[iRow];
    }
}

template <typename T>
void cuda_array_pad(const T* const arraySrc, size_t nRowSrc,
                          T* const arrayDst, size_t nRowDst,
                                             size_t nRowOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0);
    assert(nRowDst > 0);
    assert(nRowSrc + nRowOff <= nRowDst);

    cuda_array_memset(arrayDst, 0, nRowDst);

    dev_array_pad<<<iDivUp(nRowDst, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arraySrc, nRowSrc,
         arrayDst, nRowDst,
                   nRowOff);
}

// instantiation
template
void cuda_array_pad<int32_t >(const int32_t*  const arraySrc, size_t nRowSrc,
                                    int32_t*  const arrayDst, size_t nRowDst,
                                                              size_t nRowOff);
template
void cuda_array_pad<uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc,
                                    uint32_t* const arrayDst, size_t nRowDst,
                                                              size_t nRowOff);
template
void cuda_array_pad<int64_t >(const int64_t*  const arraySrc, size_t nRowSrc,
                                    int64_t*  const arrayDst, size_t nRowDst,
                                                              size_t nRowOff);
template
void cuda_array_pad<uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc,
                                    uint64_t* const arrayDst, size_t nRowDst,
                                                              size_t nRowOff);
template
void cuda_array_pad<float   >(const float*    const arraySrc, size_t nRowSrc,
                                    float*    const arrayDst, size_t nRowDst,
                                                              size_t nRowOff);
template
void cuda_array_pad<double  >(const double*   const arraySrc, size_t nRowSrc,
                                    double*   const arrayDst, size_t nRowDst,
                                                              size_t nRowOff);

// 2D
template <typename T> __global__
void dev_array_pad(const T* const arraySrc, size_t nRowSrc, size_t nColSrc,
                         T* const arrayDst, size_t nRowDst, size_t nColDst,
                                            size_t nRowOff, size_t nColOff)
{
    size_t    iRowDst = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iColDst = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowDst < nRowSrc && iColDst < nColSrc) {
        arrayDst[(iRowDst+nRowOff)*nColDst+iColDst+nColOff] = arraySrc[iRowDst*nColSrc+iColDst];
    }
}

template <typename T>
void cuda_array_pad(const T* const arraySrc, size_t nRowSrc, size_t nColSrc,
                          T* const arrayDst, size_t nRowDst, size_t nColDst,
                                             size_t nRowOff, size_t nColOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0 && nColSrc > 0);
    assert(nRowDst > 0 && nColDst > 0);
    assert(nRowSrc + nRowOff <= nRowDst);
    assert(nColSrc + nColOff <= nColDst);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nColSrc, dimBlock.x),
                     iDivUp(nRowSrc, dimBlock.y));

    cuda_array_memset(arrayDst, 0, nRowDst*nColDst);

    dev_array_pad<<<dimGrid, dimBlock>>>
        (arraySrc, nRowSrc, nColSrc,
         arrayDst, nRowDst, nColDst,
                   nRowOff, nColOff);
}
// instantiation
template
void cuda_array_pad<int32_t >(const int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                    int32_t*  const arrayDst, size_t nRowDst, size_t nColDst,
                                                              size_t nRowOff, size_t nColOff);
template
void cuda_array_pad<uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                    uint32_t* const arrayDst, size_t nRowDst, size_t nColDst,
                                                              size_t nRowOff, size_t nColOff);
template
void cuda_array_pad<int64_t >(const int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                    int64_t*  const arrayDst, size_t nRowDst, size_t nColDst,
                                                              size_t nRowOff, size_t nColOff);
template
void cuda_array_pad<uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                    uint64_t* const arrayDst, size_t nRowDst, size_t nColDst,
                                                              size_t nRowOff, size_t nColOff);
template
void cuda_array_pad<float   >(const float*    const arraySrc, size_t nRowSrc, size_t nColSrc,
                                    float*    const arrayDst, size_t nRowDst, size_t nColDst,
                                                              size_t nRowOff, size_t nColOff);
template
void cuda_array_pad<double  >(const double*   const arraySrc, size_t nRowSrc, size_t nColSrc,
                                    double*   const arrayDst, size_t nRowDst, size_t nColDst,
                                                              size_t nRowOff, size_t nColOff);

// 3D
template <typename T> __global__
void dev_array_pad(const T* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                         T* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                            size_t nRowOff, size_t nColOff, size_t nSecOff)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    size_t    iRowDst;
    size_t    iColDst = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSecDst = blockDim.x * blockIdx.x + threadIdx.x;

    if (iColDst < nColSrc && iSecDst < nSecSrc) {
        for (iRowDst = 0; iRowDst < nRowSrc; iRowDst++) {
            arrayDst[((iRowDst+nRowOff)*nColDst+(iColDst+nColOff))*nSecDst+iSecDst+nSecOff] =
            arraySrc[(iRowDst*nColSrc+iColDst)*nSecSrc+iSecDst];
        }
    }
#else
    size_t    iRowDst = blockDim.z * blockIdx.z + threadIdx.z;
    size_t    iColDst = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSecDst = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowDst < nRowSrc && iColDst < nColSrc && iSecDst < nSecSrc) {
        arrayDst[((iRowDst+nRowOff)*nColDst+(iColDst+nColOff))*nSecDst+iSecDst+nSecOff] =
        arraySrc[(iRowDst*nColSrc+iColDst)*nSecSrc+iSecDst];
    }
#endif
}

template <typename T>
void cuda_array_pad(const T* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                          T* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                             size_t nRowOff, size_t nColOff, size_t nSecOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);
    assert(nRowDst > 0 && nColDst > 0 && nSecDst > 0);
    assert(nRowSrc + nRowOff <= nRowDst);
    assert(nColSrc + nColOff <= nColDst);
    assert(nSecSrc + nSecOff <= nSecDst);

#ifdef __GEM_CUDA_ARCH_HOST_130__
    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nSecSrc, dimBlock.x),
                     iDivUp(nColSrc, dimBlock.y));
#else
    dim3    dimBlock(BLOCK_3D_NSEC, BLOCK_3D_NCOL, BLOCK_3D_NROW);
    dim3    dimGrid (iDivUp(nSecSrc, dimBlock.x),
                     iDivUp(nColSrc, dimBlock.y),
                     iDivUp(nRowSrc, dimBlock.z));
#endif

    cuda_array_memset(arrayDst, 0, nRowDst*nColDst*nSecDst);

    dev_array_pad<<<dimGrid, dimBlock>>>
        (arraySrc, nRowSrc, nColSrc, nSecSrc,
         arrayDst, nRowDst, nColDst, nSecDst,
                    nRowOff, nColOff, nSecOff);
}

// instantiation
template
void cuda_array_pad<int32_t >(const int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                    int32_t*  const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                              size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_pad<uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                    uint32_t* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                              size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_pad<int64_t >(const int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                    int64_t*  const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                              size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_pad<uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                    uint64_t* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                              size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_pad<float   >(const float*    const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                    float*    const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                              size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_pad<double  >(const double*   const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                    double*   const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                              size_t nRowOff, size_t nColOff, size_t nSecOff);

} // namespace gem

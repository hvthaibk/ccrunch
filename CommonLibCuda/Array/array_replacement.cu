/***********************************************************************
 *  File:       array_replacement.cu
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
 * Replacement
 ****************************************/

// 1D
template <typename T1, typename T2> __global__
void dev_array_replace(      T1* const arraySrc, size_t nRowSrc,
                       const T2* const arrayRep, size_t nRowRep,
                                                 size_t nRowOff)
{
    size_t    iRowRep = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowRep < nRowRep) {
        arraySrc[iRowRep+nRowOff] = (T1) arrayRep[iRowRep];
    }
}

template <typename T1, typename T2>
void cuda_array_replace(      T1* const arraySrc, size_t nRowSrc,
                        const T2* const arrayRep, size_t nRowRep,
                                                  size_t nRowOff)
{
    assert(arraySrc != NULL);
    assert(arrayRep != NULL);
    assert(nRowSrc > 0);
    assert(nRowRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);

    dev_array_replace<<<iDivUp(nRowRep, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arraySrc, nRowSrc,
         arrayRep, nRowRep,
                   nRowOff);
}

// instantiation
template
void cuda_array_replace<int32_t ,int32_t >(      int32_t*  const arraySrc, size_t nRowSrc,
                                           const int32_t*  const arrayRep, size_t nRowRep,
                                                                           size_t nRowOff);
template
void cuda_array_replace<uint32_t,uint32_t>(      uint32_t* const arraySrc, size_t nRowSrc,
                                           const uint32_t* const arrayRep, size_t nRowRep,
                                                                           size_t nRowOff);
template
void cuda_array_replace<int64_t ,int64_t >(      int64_t*  const arraySrc, size_t nRowSrc,
                                           const int64_t*  const arrayRep, size_t nRowRep,
                                                                           size_t nRowOff);
template
void cuda_array_replace<uint64_t,uint64_t>(      uint64_t* const arraySrc, size_t nRowSrc,
                                           const uint64_t* const arrayRep, size_t nRowRep,
                                                                           size_t nRowOff);
template
void cuda_array_replace<float   ,float   >(      float*    const arraySrc, size_t nRowSrc,
                                           const float*    const arrayRep, size_t nRowRep,
                                                                           size_t nRowOff);
template
void cuda_array_replace<double  ,double  >(      double*   const arraySrc, size_t nRowSrc,
                                           const double*   const arrayRep, size_t nRowRep,
                                                                           size_t nRowOff);

// 1D
template <typename T1, typename T2> __global__
void dev_array_replace(T1* const arraySrc, size_t nRowSrc,
                       T2        value,    size_t nRowRep,
                                           size_t nRowOff)
{
    size_t    iRowRep = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowRep < nRowRep) {
        arraySrc[iRowRep+nRowOff] = (T1) value;
    }
}

template <typename T1, typename T2>
void cuda_array_replace(T1* const arraySrc, size_t nRowSrc,
                        T2        value,    size_t nRowRep,
                                            size_t nRowOff)
{
    assert(arraySrc != NULL);
    assert(nRowSrc > 0);
    assert(nRowRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);

    dev_array_replace<<<iDivUp(nRowRep, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arraySrc, nRowSrc,
         value,    nRowRep,
                   nRowOff);
}

// instantiation
template
void cuda_array_replace<int32_t ,int32_t >(int32_t*  const arraySrc, size_t nRowSrc,
                                           int32_t         value,    size_t nRowRep,
                                                                     size_t nRowOff);
template
void cuda_array_replace<uint32_t,uint32_t>(uint32_t* const arraySrc, size_t nRowSrc,
                                           uint32_t        value,    size_t nRowRep,
                                                                     size_t nRowOff);
template
void cuda_array_replace<int64_t ,int64_t >(int64_t*  const arraySrc, size_t nRowSrc,
                                           int64_t         value,    size_t nRowRep,
                                                                     size_t nRowOff);
template
void cuda_array_replace<uint64_t,uint64_t>(uint64_t* const arraySrc, size_t nRowSrc,
                                           uint64_t        value,    size_t nRowRep,
                                                                     size_t nRowOff);
template
void cuda_array_replace<float   ,float   >(float*    const arraySrc, size_t nRowSrc,
                                           float           value,    size_t nRowRep,
                                                                     size_t nRowOff);
template
void cuda_array_replace<double  ,double  >(double*   const arraySrc, size_t nRowSrc,
                                           double          value,    size_t nRowRep,
                                                                     size_t nRowOff);

// 2D
template <typename T1, typename T2> __global__
void dev_array_replace(      T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                       const T2* const arrayRep, size_t nRowRep, size_t nColRep,
                                                 size_t nRowOff, size_t nColOff)
{
    size_t    iRowRep = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iColRep = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowRep < nRowRep && iColRep < nColRep) {
        arraySrc[(iRowRep+nRowOff)*nColSrc+iColRep+nColOff] =
        (T1) arrayRep[iRowRep*nColRep+iColRep];
    }
}

template <typename T1, typename T2>
void cuda_array_replace(      T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                        const T2* const arrayRep, size_t nRowRep, size_t nColRep,
                                                  size_t nRowOff, size_t nColOff)
{
    assert(arraySrc != NULL);
    assert(arrayRep != NULL);
    assert(nRowSrc > 0 && nColSrc > 0);
    assert(nRowRep > 0 && nColRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);
    assert(nColRep + nColOff <= nColSrc);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nColRep, dimBlock.x),
                     iDivUp(nRowRep, dimBlock.y));

    dev_array_replace<<<dimGrid, dimBlock>>>
        (arraySrc, nRowSrc, nColSrc,
         arrayRep, nRowRep, nColRep,
                   nRowOff, nColOff);
}

// instantiation
template
void cuda_array_replace<int32_t ,int32_t >(      int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           const int32_t*  const arrayRep, size_t nRowRep, size_t nColRep,
                                                                           size_t nRowOff, size_t nColOff);
template
void cuda_array_replace<uint32_t,uint32_t>(      uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           const uint32_t* const arrayRep, size_t nRowRep, size_t nColRep,
                                                                           size_t nRowOff, size_t nColOff);
template
void cuda_array_replace<int64_t ,int64_t >(      int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           const int64_t*  const arrayRep, size_t nRowRep, size_t nColRep,
                                                                           size_t nRowOff, size_t nColOff);
template
void cuda_array_replace<uint64_t,uint64_t>(      uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           const uint64_t* const arrayRep, size_t nRowRep, size_t nColRep,
                                                                           size_t nRowOff, size_t nColOff);
template
void cuda_array_replace<float   ,float   >(      float*    const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           const float*    const arrayRep, size_t nRowRep, size_t nColRep,
                                                                           size_t nRowOff, size_t nColOff);
template
void cuda_array_replace<double  ,double  >(      double*   const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           const double*   const arrayRep, size_t nRowRep, size_t nColRep,
                                                                           size_t nRowOff, size_t nColOff);

// 2D
template <typename T1, typename T2> __global__
void dev_array_replace(T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                       T2        value,    size_t nRowRep, size_t nColRep,
                                           size_t nRowOff, size_t nColOff)
{
    size_t    iRowRep = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iColRep = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowRep < nRowRep && iColRep < nColRep) {
        arraySrc[(iRowRep+nRowOff)*nColSrc+iColRep+nColOff] = (T1) value;
    }
}

template <typename T1, typename T2>
void cuda_array_replace(T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                        T2        value,    size_t nRowRep, size_t nColRep,
                                            size_t nRowOff, size_t nColOff)
{
    assert(arraySrc != NULL);
    assert(nRowSrc > 0 && nColSrc > 0);
    assert(nRowRep > 0 && nColRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);
    assert(nColRep + nColOff <= nColSrc);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nColRep, dimBlock.x),
                     iDivUp(nRowRep, dimBlock.y));

    dev_array_replace<<<dimGrid, dimBlock>>>
        (arraySrc, nRowSrc, nColSrc,
         value,    nRowRep, nColRep,
                   nRowOff, nColOff);
}

// instantiation
template
void cuda_array_replace<int32_t ,int32_t >(int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           int32_t         value,    size_t nRowRep, size_t nColRep,
                                                                     size_t nRowOff, size_t nColOff);
template
void cuda_array_replace<uint32_t,uint32_t>(uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           uint32_t        value,    size_t nRowRep, size_t nColRep,
                                                                     size_t nRowOff, size_t nColOff);
template
void cuda_array_replace<int64_t ,int64_t >(int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           int64_t         value,    size_t nRowRep, size_t nColRep,
                                                                     size_t nRowOff, size_t nColOff);
template
void cuda_array_replace<uint64_t,uint64_t>(uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           uint64_t        value,    size_t nRowRep, size_t nColRep,
                                                                     size_t nRowOff, size_t nColOff);
template
void cuda_array_replace<float   ,float   >(float*    const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           float           value,    size_t nRowRep, size_t nColRep,
                                                                     size_t nRowOff, size_t nColOff);
template
void cuda_array_replace<double  ,double  >(double*   const arraySrc, size_t nRowSrc, size_t nColSrc,
                                           double          value,    size_t nRowRep, size_t nColRep,
                                                                     size_t nRowOff, size_t nColOff);

// 3D
template <typename T1, typename T2> __global__
void dev_array_replace(      T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                       const T2* const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                 size_t nRowOff, size_t nColOff, size_t nSecOff)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    size_t    iRowRep;
    size_t    iColRep = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSecRep = blockDim.x * blockIdx.x + threadIdx.x;

    if (iColRep < nColRep && iSecRep < nSecRep) {
        for (iRowRep = 0; iRowRep < nRowRep; iRowRep++) {
            arraySrc[((iRowRep+nRowOff)*nColSrc+iColRep+nColOff)*nSecSrc+iSecRep+nSecOff] =
            (T1) arrayRep[(iRowRep*nColRep+iColRep)*nSecRep+iSecRep];
        }
    }
#else
    size_t    iRowRep = blockDim.z * blockIdx.z + threadIdx.z;
    size_t    iColRep = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSecRep = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowRep < nRowRep && iColRep < nColRep && iSecRep < nSecRep) {
        arraySrc[((iRowRep+nRowOff)*nColSrc+iColRep+nColOff)*nSecSrc+iSecRep+nSecOff] =
        (T1) arrayRep[(iRowRep*nColRep+iColRep)*nSecRep+iSecRep];
    }
#endif
}

template <typename T1, typename T2>
void cuda_array_replace(      T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                        const T2* const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                  size_t nRowOff, size_t nColOff, size_t nSecOff)
{
    assert(arraySrc != NULL);
    assert(arrayRep != NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);
    assert(nRowRep > 0 && nColRep > 0 && nSecRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);
    assert(nColRep + nColOff <= nColSrc);
    assert(nSecRep + nSecOff <= nSecSrc);

#ifdef __GEM_CUDA_ARCH_HOST_130__
    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nSecRep, dimBlock.x),
                     iDivUp(nColRep, dimBlock.y));
#else
    dim3    dimBlock(BLOCK_3D_NSEC, BLOCK_3D_NCOL, BLOCK_3D_NROW);
    dim3    dimGrid (iDivUp(nSecRep, dimBlock.x),
                     iDivUp(nColRep, dimBlock.y),
                     iDivUp(nRowRep, dimBlock.z));
#endif

    dev_array_replace<<<dimGrid, dimBlock>>>
        (arraySrc, nRowSrc, nColSrc, nSecSrc,
         arrayRep, nRowRep, nColRep, nSecRep,
                   nRowOff, nColOff, nSecOff);
}

// instantiation
template
void cuda_array_replace<int32_t ,int32_t >(      int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           const int32_t*  const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                           size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_replace<uint32_t,uint32_t>(      uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           const uint32_t* const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                           size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_replace<int64_t ,int64_t >(      int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           const int64_t*  const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                           size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_replace<uint64_t,uint64_t>(      uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           const uint64_t* const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                           size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_replace<float   ,float   >(      float*    const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           const float*    const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                           size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_replace<double  ,double  >(      double*   const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           const double*   const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                           size_t nRowOff, size_t nColOff, size_t nSecOff);

// 3D
template <typename T1, typename T2> __global__
void dev_array_replace(T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                       T2        value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                           size_t nRowOff, size_t nColOff, size_t nSecOff)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    size_t    iRowRep;
    size_t    iColRep = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSecRep = blockDim.x * blockIdx.x + threadIdx.x;

    if (iColRep < nColRep && iSecRep < nSecRep) {
        for (iRowRep = 0; iRowRep < nRowRep; iRowRep++) {
            arraySrc[((iRowRep+nRowOff)*nColSrc+iColRep+nColOff)*nSecSrc+iSecRep+nSecOff] = (T1) value;
        }
    }
#else
    size_t    iRowRep = blockDim.z * blockIdx.z + threadIdx.z;
    size_t    iColRep = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSecRep = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowRep < nRowRep && iColRep < nColRep && iSecRep < nSecRep) {
        arraySrc[((iRowRep+nRowOff)*nColSrc+iColRep+nColOff)*nSecSrc+iSecRep+nSecOff] = (T1) value;
    }
#endif
}

template <typename T1, typename T2>
void cuda_array_replace(T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                        T2        value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                            size_t nRowOff, size_t nColOff, size_t nSecOff)
{
    assert(arraySrc != NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);
    assert(nRowRep > 0 && nColRep > 0 && nSecRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);
    assert(nColRep + nColOff <= nColSrc);
    assert(nSecRep + nSecOff <= nSecSrc);

#ifdef __GEM_CUDA_ARCH_HOST_130__
    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nSecRep, dimBlock.x),
                     iDivUp(nColRep, dimBlock.y));
#else
    dim3    dimBlock(BLOCK_3D_NSEC, BLOCK_3D_NCOL, BLOCK_3D_NROW);
    dim3    dimGrid (iDivUp(nSecRep, dimBlock.x),
                     iDivUp(nColRep, dimBlock.y),
                     iDivUp(nRowRep, dimBlock.z));
#endif

    dev_array_replace<<<dimGrid, dimBlock>>>
        (arraySrc, nRowSrc, nColSrc, nSecSrc,
         value,    nRowRep, nColRep, nSecRep,
                   nRowOff, nColOff, nSecOff);
}

// instantiation
template
void cuda_array_replace<int32_t ,int32_t >(int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           int32_t         value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                     size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_replace<uint32_t,uint32_t>(uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           uint32_t        value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                     size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_replace<int64_t ,int64_t >(int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           int64_t         value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                     size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_replace<uint64_t,uint64_t>(uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           uint64_t        value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                     size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_replace<float   ,float   >(float*    const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           float           value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                     size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void cuda_array_replace<double  ,double  >(double*   const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                           double          value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                     size_t nRowOff, size_t nColOff, size_t nSecOff);

} // namespace gem

/***********************************************************************
 *  File:       array.cuh
 *
 *  Purpose:    Header file for array-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_DEV_ARRAY_CUH__
#define __GEM_DEV_ARRAY_CUH__

#include "array.hpp"

#include <cuda.h>
#include <cuComplex.h>
#include <cuda_runtime.h>

namespace gem {

const unsigned int BLOCK_1D_NROW = 512;
const unsigned int BLOCK_2D_NROW = 16;
const unsigned int BLOCK_2D_NCOL = 16;
const unsigned int BLOCK_3D_NROW = 8;
const unsigned int BLOCK_3D_NCOL = 8;
const unsigned int BLOCK_3D_NSEC = 8;

const unsigned int TILE_DIM   = 32;
const unsigned int BLOCK_ROWS = 8;

// assert() is only supported for devices of compute capability >= 2.0
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)) || defined(__GEM_CUDA_NVCC_VERSION_40__)
#undef  assert
#define assert(arg)
#endif

template<typename T>
struct SharedMemory
{
    __device__ inline operator       T*()
    {
        extern __shared__ int32_t __smem[];
        return (T*)__smem;
    }

    __device__ inline operator const T*() const
    {
        extern __shared__ int32_t __smem[];
        return (T*)__smem;
    }
};

template<>
struct SharedMemory<int64_t>
{
    __device__ inline operator       int64_t*()
    {
        extern __shared__ int64_t __smem_d[];
        return (int64_t*)__smem_d;
    }

    __device__ inline operator const int64_t*() const
    {
        extern __shared__ int64_t __smem_d[];
        return (int64_t*)__smem_d;
    }
};

template<>
struct SharedMemory<uint64_t>
{
    __device__ inline operator       uint64_t*()
    {
        extern __shared__ int64_t __smem_d[];
        return (uint64_t*)__smem_d;
    }

    __device__ inline operator const uint64_t*() const
    {
        extern __shared__ int64_t __smem_d[];
        return (uint64_t*)__smem_d;
    }
};

template<>
struct SharedMemory<double>
{
    __device__ inline operator       double*()
    {
        extern __shared__ int64_t __smem_d[];
        return (double*)__smem_d;
    }

    __device__ inline operator const double*() const
    {
        extern __shared__ int64_t __smem_d[];
        return (double*)__smem_d;
    }
};

/*****************************************
 * memcpy / memset / memcheck
 ****************************************/

// memcpy
template <typename T>
void cuda_array_memcpy_d2d(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void cuda_array_memcpy_d2h(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void cuda_array_memcpy_h2d(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void cuda_array_memcpy_h2h(T* const arrayDst, const T* const arraySrc, size_t nRow);

void cuda_array_memcpy_d2d(cudaArray* const arrayDst, const float* const arraySrc, size_t nRow, size_t nCol);
void cuda_array_memcpy_h2d(cudaArray* const arrayDst, const float* const arraySrc, size_t nRow, size_t nCol);

void cuda_array_memcpy_d2d(cudaArray* const arrayDst, const float* const arraySrc, size_t nRow, size_t nCol, size_t nSec);
void cuda_array_memcpy_h2d(cudaArray* const arrayDst, const float* const arraySrc, size_t nRow, size_t nCol, size_t nSec);

// memset
template <typename T>
void cuda_array_memset(T* const array, int value, size_t nRow);

// check memory leak
void cuda_array_memcheck(void);

/*****************************************
 * Array allocation/deallocation
 *       w/o zero initialization
 ****************************************/

// cudaArray
void cuda_arrayArr_new   (cudaArray*& array, size_t nRow, size_t nCol, size_t nSec);
void cuda_arrayArr_new   (cudaArray*& array, size_t nRow, size_t nCol);
void cuda_arrayArr_delete(cudaArray*& array);

// device memory
// size(array1D) = [nRow]
template <typename T>
void cuda_arrayDev_new(T*& array1D, size_t nRow);

template <typename T>
void cuda_arrayDev_new_zero(T*& array1D, size_t nRow);

template <typename T>
void cuda_arrayDev_delete(T*& array1D);

// size(array2D) = [nRow nCol]
template <typename T>
void cuda_arrayDev_new(T**& array2D, size_t nRow, size_t nCol);

template <typename T>
void cuda_arrayDev_new_zero(T**& array2D, size_t nRow, size_t nCol);

template <typename T>
void cuda_arrayDev_delete(T**& array2D, size_t nRow);

// host memory
// size(array1D) = [nRow]
template <typename T>
void cuda_arrayHost_new(T*& array1D, size_t nRow);

template <typename T>
void cuda_arrayHost_new_zero(T*& array1D, size_t nRow);

template <typename T>
void cuda_arrayHost_delete(T*& array1D);

/*****************************************
 * Padding
 ****************************************/

// 1D
template <typename T>
void cuda_array_pad(const T* const arraySrc, size_t nRowSrc,
                          T* const arrayDst, size_t nRowDst,
                                             size_t nRowOff);

// 2D
template <typename T>
void cuda_array_pad(const T* const arraySrc, size_t nRowSrc, size_t nColSrc,
                          T* const arrayDst, size_t nRowDst, size_t nColDst,
                                             size_t nRowOff, size_t nColOff);

// 3D
template <typename T>
void cuda_array_pad(const T* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                          T* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                             size_t nRowOff, size_t nColOff, size_t nSecOff);


/*****************************************
 * Cropping
 ****************************************/

// 1D
template <typename T>
void cuda_array_crop(const T* const arraySrc, size_t nRowSrc,
                           T* const arrayDst, size_t nRowDst,
                                              size_t nRowOff);

// 2D
template <typename T>
void cuda_array_crop(const T* const arraySrc, size_t nRowSrc, size_t nColSrc,
                           T* const arrayDst, size_t nRowDst, size_t nColDst,
                                              size_t nRowOff, size_t nColOff);

// 3D
template <typename T>
void cuda_array_crop(const T* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                           T* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                              size_t nRowOff, size_t nColOff, size_t nSecOff);

/*****************************************
 * Replacement
 ****************************************/

// 1D
template <typename T1, typename T2>
void cuda_array_replace(T1* const arraySrc, size_t nRowSrc,
                        T2        value,    size_t nRowRep,
                                            size_t nRowOff);
// 2D
template <typename T1, typename T2>
void cuda_array_replace(T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                        T2        value,    size_t nRowRep, size_t nColRep,
                                            size_t nRowOff, size_t nColOff);
// 3D
template <typename T1, typename T2>
void cuda_array_replace(T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                        T2        value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                            size_t nRowOff, size_t nColOff, size_t nSecOff);

// 1D
template <typename T1, typename T2>
void cuda_array_replace(      T1* const arraySrc, size_t nRowSrc,
                        const T2* const arrayRep, size_t nRowRep,
                                                  size_t nRowOff);
// 2D
template <typename T1, typename T2>
void cuda_array_replace(      T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                        const T2* const arrayRep, size_t nRowRep, size_t nColRep,
                                                  size_t nRowOff, size_t nColOff);
// 3D
template <typename T1, typename T2>
void cuda_array_replace(      T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                        const T2* const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                  size_t nRowOff, size_t nColOff, size_t nSecOff);

/*****************************************
 * Indexing
 ****************************************/

// 1D index -> 2D sub-indices
template <typename T>
void cuda_array_index_ind2sub(const T* const arrayInd,
                                    T* const arrayRow,
                                    T* const arrayCol,
                                    size_t   length,
                                    size_t   nRow, size_t nCol);

// 1D index -> 3D sub-indices
template <typename T>
void cuda_array_index_ind2sub(const T* const arrayInd,
                                    T* const arrayRow,
                                    T* const arrayCol,
                                    T* const arraySec,
                                    size_t   length,
                                    size_t   nRow, size_t nCol, size_t nSec);

// 2D sub-indices -> 1D index
template <typename T>
void cuda_array_index_sub2ind(const T* const arrayRow,
                              const T* const arrayCol,
                                    T* const arrayInd,
                                    size_t   length,
                                    size_t   nRow, size_t nCol);

// 3D sub-indices -> 1D index
template <typename T>
void cuda_array_index_sub2ind(const T* const arrayRow,
                              const T* const arrayCol,
                              const T* const arraySec,
                                    T* const arrayInd,
                                    size_t   length,
                                    size_t   nRow, size_t nCol, size_t nSec);

// column-major -> row-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow)
template <typename T1, typename T2>
void cuda_array_index_col2row(const T1* const arrayCol, T2* const arrayRow,
                              size_t nRow, size_t nCol);

// column-major -> row-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow), nSec = dim3(arrayRow)
template <typename T1, typename T2>
void cuda_array_index_col2row(const T1* const arrayCol, T2* const arrayRow,
                              size_t nRow, size_t nCol, size_t nSec);

// row-major -> column-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow)
template <typename T1, typename T2>
void cuda_array_index_row2col(const T1* const arrayRow, T2* const arrayCol,
                              size_t nRow, size_t nCol);

// row-major -> column-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow), nSec = dim3(arrayRow)
template <typename T1, typename T2>
void cuda_array_index_row2col(const T1* const arrayRow, T2* const arrayCol,
                              size_t nRow, size_t nCol, size_t nSec);

/*****************************************
 * Permutation
 ****************************************/

// 2D
template <typename T>
void cuda_array_permute(const T* const arraySrc, T* const arrayDst,
                        size_t nRow, size_t nCol);

// 3D
template <typename T>
void cuda_array_permute(const T* const arraySrc, T* const arrayDst,
                        size_t nRow, size_t nCol, size_t nSec,
                        ePermute permute);

/*****************************************
 * Element-wise operations
 ****************************************/

// ismember
template <typename T>
void cuda_array_ismember(const T* const arraySrc, size_t nRowSrc,
                         const T* const arraySet, size_t nRowSet,
                               T* const arrayMem);

// set value
template <typename T>
void cuda_array_setval(T* const array, T value, size_t nRow);

// type cast
template <typename T1, typename T2>
void cuda_array_typecast(const T1* const arraySrc, size_t nRow,
                               T2* const arrayDst);

// threshold elements to 0/1
template <typename T>
void cuda_array_threshold(T* const array, size_t nRow,
                          T        cutoff);
template <typename T1, typename T2>
void cuda_array_threshold(const T1* const arraySrc, size_t nRow,
                                T1        cutoff,
                                T2* const arrayDst);
template <typename T1, typename T2>
void cuda_array_threshold(const T1* const arraySrc, size_t nRow,
                          const T1* const arrayThr,
                                T2* const arrayDst);

// cutoff small elements to value
template <typename T>
void cuda_array_cutoff(T* const array, size_t nRow,
                       T        cutoff,
                       T        value);
template <typename T>
void cuda_array_cutoff(const T* const arraySrc, size_t nRow,
                             T        cutoff,
                             T* const arrayDst,
                             T        value);
template <typename T>
void cuda_array_cutoff(const T* const arraySrc, size_t nRow,
                       const T* const arrayThr,
                             T* const arrayDst,
                             T        value);

// suppress big elements to value
template <typename T>
void cuda_array_suppress(T* const array, size_t nRow,
                         T        cutoff,
                         T        value);
template <typename T>
void cuda_array_suppress(const T* const arraySrc, size_t nRow,
                               T        cutoff,
                               T* const arrayDst,
                               T        value);
template <typename T>
void cuda_array_suppress(const T* const arraySrc, size_t nRow,
                         const T* const arrayThr,
                               T* const arrayDst,
                               T        value);

// remove zero elements by adding noise
template <typename T>
void cuda_array_remove_rezo(const T* const arrayNoise, size_t nRow,
                                  T* const arrayData,
                                  T        cutoff);

/*****************************************
 * Math
 ****************************************/

// abs
template <typename T>
void cuda_array_math_abs(T* const array, size_t nRow);
template <typename T>
void cuda_array_math_abs(T* const arrayDst, const T* const arraySrc, size_t nRow);

// inverse
template <typename T>
void cuda_array_math_inv(T* const array, size_t nRow);
template <typename T>
void cuda_array_math_inv(T* const arrayDst, const T* const arraySrc, size_t nRow);

// norm (mean = 0, std = 1)
template <typename T>
void cuda_array_math_norm(T* const array, size_t nRow);
template <typename T>
void cuda_array_math_norm(T* const arrayDst, const T* const arraySrc, size_t nRow);

// square
template <typename T>
void cuda_array_math_sqr(T* const array, size_t nRow);
template <typename T>
void cuda_array_math_sqr(T* const arrayDst, const T* const arraySrc, size_t nRow);

// square root
template <typename T>
void cuda_array_math_sqrt(T* const array, size_t nRow);
template <typename T>
void cuda_array_math_sqrt(T* const arrayDst, const T* const arraySrc, size_t nRow);

// add
template <typename T>
void cuda_array_math_add(T* const array, T value, size_t nRow);
template <typename T>
void cuda_array_math_add(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void cuda_array_math_add(T* const arrayDst, const T* const arraySrc, T value, size_t nRow);
template <typename T>
void cuda_array_math_add(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow);

// substract
template <typename T>
void cuda_array_math_sub(T* const array, T value, size_t nRow);
template <typename T>
void cuda_array_math_sub(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void cuda_array_math_sub(T* const arrayDst, const T* const arraySrc, T value, size_t nRow);
template <typename T>
void cuda_array_math_sub(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow);

// multiply
template <typename T>
void cuda_array_math_mul(T* const array, T value, size_t nRow);
template <typename T>
void cuda_array_math_mul(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void cuda_array_math_mul(T* const arrayDst, const T* const arraySrc, T value, size_t nRow);
template <typename T>
void cuda_array_math_mul(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow);

// divide
template <typename T>
void cuda_array_math_div(T* const array, T value, size_t nRow);
template <typename T>
void cuda_array_math_div(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void cuda_array_math_div(T* const arrayDst, const T* const arraySrc, T value, size_t nRow);
template <typename T>
void cuda_array_math_div(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow);

/*****************************************
 * BLAS
 ****************************************/

// blas axpy: Y[i] = alpha*X[i] + Y[i]
template <typename T>
void cuda_array_blas_axpy(      T* const arrayY,
                          const T* const arrayX,
                          T alpha, size_t nRow);

// blas gemm: C = A*B, nColA = nRowB
template <typename T>
void cuda_array_blas_gemm(      T* const arrayC,
                          const T* const arrayA,
                          const T* const arrayB,
                          size_t nRowA, size_t nColA, size_t nColB,
                          unsigned int mode = 0);

/*****************************************
 * Reduction
 ****************************************/

// max, min, sum, sum of squares
template <typename T>
T cuda_array_reduce_max(const T* const array, size_t nRow, unsigned int nblocksCUDA = 256);
template <typename T>
T cuda_array_reduce_min(const T* const array, size_t nRow, unsigned int nblocksCUDA = 256);
template <typename T>
T cuda_array_reduce_sum(const T* const array, size_t nRow, unsigned int nblocksCUDA = 256);
template <typename T>
T cuda_array_reduce_sum2(const T* const array, size_t nRow, unsigned int nblocksCUDA = 256);

// mean, std
template <typename T>
T cuda_array_reduce_mean(const T* const array, size_t nRow, unsigned int nblocksCUDA = 256);
template <typename T>
T cuda_array_reduce_std(const T* const array, size_t nRow, unsigned int nblocksCUDA = 256);

// index
template <typename T>
size_t cuda_array_reduce_maxindx(const T* const array, size_t nRow, unsigned int nblocksCUDA = 256);
template <typename T>
size_t cuda_array_reduce_minindx(const T* const array, size_t nRow, unsigned int nblocksCUDA = 256);

// corr, conv
template <typename T>
T cuda_array_reduce_corr(const T* const arraySrc1,
                         const T* const arraySrc2,
                         size_t nRow, unsigned int nblocksCUDA = 256);
template <typename T>
T cuda_array_reduce_conv(const T* const arraySrc1,
                         const T* const arraySrc2,
                         size_t nRow, unsigned int nblocksCUDA = 256);

// thrust
template <typename T>
T cuda_array_reduce_thrust(const T* const array, size_t nRow, eReduce reduce);

/*****************************************
 * Reduction NDim
 ****************************************/

template <typename T>
void cuda_array_reducendim_max (const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_ROW);
template <typename T>
void cuda_array_reducendim_min (const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_ROW);
template <typename T>
void cuda_array_reducendim_sum (const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_ROW);
template <typename T>
void cuda_array_reducendim_sum2(const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_ROW);
template <typename T>
void cuda_array_reducendim_mean(const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_ROW);
template <typename T>
void cuda_array_reducendim_std (const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_ROW);

template <typename T>
void cuda_array_reducendim_max (const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_SEC);
template <typename T>
void cuda_array_reducendim_min (const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_SEC);
template <typename T>
void cuda_array_reducendim_sum (const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_SEC);
template <typename T>
void cuda_array_reducendim_sum2(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_SEC);
template <typename T>
void cuda_array_reducendim_mean(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_SEC);
template <typename T>
void cuda_array_reducendim_std(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim = REDUCE_NDIM_SEC);

/*****************************************
 * Circular shifting
 ****************************************/

// 1D
template <typename T>
void cuda_array_circshift(const T* const arraySrc, T* const arrayDst,
                          size_t    nRow,
                          ptrdiff_t nRowOff);

// 2D
template <typename T>
void cuda_array_circshift(const T* const arraySrc, T* const arrayDst,
                          size_t    nRow,    size_t    nCol,
                          ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D
template <typename T>
void cuda_array_circshift(const T* const arraySrc, T* const arrayDst,
                          size_t    nRow,    size_t    nCol,    size_t    nSec,
                          ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);


/*****************************************
 * Shifting
 ****************************************/

// 1D
template <typename T>
void cuda_array_shift(const T* const arraySrc, T* const arrayDst,
                      size_t    nRow,
                      ptrdiff_t nRowOff);

// 2D
template <typename T>
void cuda_array_shift(const T* const arraySrc, T* const arrayDst,
                      size_t    nRow,    size_t    nCol,
                      ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D
template <typename T>
void cuda_array_shift(const T* const arraySrc, T* const arrayDst,
                      size_t    nRow,    size_t    nCol,    size_t    nSec,
                      ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

/*****************************************
 * Copying
 ****************************************/

// 1D
template <typename T>
void cuda_array_copy(const T* const arraySrc, T* const arrayDst, size_t nRow);

// 2D
template <typename T>
void cuda_array_copy(const T* const arraySrc, T* const arrayDst, size_t nRow, size_t nCol);

// 3D
template <typename T>
void cuda_array_copy(const T* const arraySrc, T* const arrayDst, size_t nRow, size_t nCol, size_t nSec);

/*****************************************
 * Mask
 ****************************************/

// 1D copy
template <typename T>
void cuda_array_mask_copy(const T* const array,  size_t    nRow,
                          const T* const mask,   size_t    nRowMsk,
                                T* const output, ptrdiff_t nRowOff);

// 2D copy
template <typename T>
void cuda_array_mask_copy(const T* const array,  size_t    nRow,    size_t    nCol,
                          const T* const mask,   size_t    nRowMsk, size_t    nColMsk,
                                T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D copy
template <typename T>
void cuda_array_mask_copy(const T* const array,  size_t    nRow,    size_t    nCol,    size_t    nSec,
                          const T* const mask,   size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

// 1D replace
template <typename T>
void cuda_array_mask_replace(      T* const array, size_t    nRow,
                             const T* const mask,  size_t    nRowMsk,
                                   T        value, ptrdiff_t nRowOff,
                                   bool     maskInside = true);

// 2D replace
template <typename T>
void cuda_array_mask_replace(      T* const array, size_t    nRow,    size_t    nCol,
                             const T* const mask,  size_t    nRowMsk, size_t    nColMsk,
                                   T        value, ptrdiff_t nRowOff, ptrdiff_t nColOff,
                                   bool     maskInside = true);
// 3D replace
template <typename T>
void cuda_array_mask_replace(      T* const array, size_t    nRow,    size_t    nCol,    size_t    nSec,
                             const T* const mask,  size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                   T        value, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff,
                                   bool     maskInside = true);

// sum
template <typename T>
T cuda_array_mask_sum(const T* const array, const T* const mask, size_t nRow, unsigned int nblocksCUDA = 256);

// substract
template <typename T>
void cuda_array_mask_sub(T* const arrayDst, const T* const arraySrc, const T* const mask, T value, size_t nRow);

// substract - square - sum
template <typename T>
T cuda_array_mask_subsqrsum(const T* const array, T value, const T* const mask, size_t nRow, unsigned int nblocksCUDA = 256);

// substract - divide - multiply
template <typename T>
void cuda_array_mask_subdivmul(T* const array, T valsub, T valdiv, const T* const mask, size_t nRow);

/*****************************************
 * Random
 ****************************************/

// normal distribution
template <typename T>
void cuda_array_random_normal(T* const array, size_t nRow, T mean, T std, size_t seed = 32);

// uniform distribution
template <typename T>
void cuda_array_random_uniform(T* const array, size_t nRow, T min, T max, size_t seed = 32);

} // namespace gem

#endif

/***********************************************************************
 *  File:       array_math.cu
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
 * Math
 ****************************************/

// abs
template <typename T> __global__
void dev_array_math_abs(T* const array, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] = fabsf(array[iRow]);
    }
}

template <> __global__
void dev_array_math_abs<int32_t>(int32_t* const array, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] = abs(array[iRow]);
    }
}

template <> __global__
void dev_array_math_abs<uint32_t>(uint32_t* const array, size_t nRow) { }

template <> __global__
void dev_array_math_abs<int64_t>(int64_t* const array, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] = abs(array[iRow]);
    }
}

template <> __global__
void dev_array_math_abs<uint64_t>(uint64_t* const array, size_t nRow) { }

template <> __global__
void dev_array_math_abs<double>(double* const array, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] = fabs(array[iRow]);
    }
}

template <typename T> __global__
void dev_array_math_abs(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = fabsf(arraySrc[iRow]);
    }
}

template <> __global__
void dev_array_math_abs<int32_t>(int32_t* const arrayDst, const int32_t* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = abs(arraySrc[iRow]);
    }
}

template <> __global__
void dev_array_math_abs<uint32_t>(uint32_t* const , const uint32_t* const , size_t ) { }

template <> __global__
void dev_array_math_abs<int64_t>(int64_t* const arrayDst, const int64_t* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = abs(arraySrc[iRow]);
    }
}

template <> __global__
void dev_array_math_abs<uint64_t>(uint64_t* const , const uint64_t* const , size_t ) { }

template <> __global__
void dev_array_math_abs<double>(double* const arrayDst, const double* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = fabs(arraySrc[iRow]);
    }
}

template <typename T>
void cuda_array_math_abs(T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    dev_array_math_abs
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array, nRow);
}

template <typename T>
void cuda_array_math_abs(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_abs
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, nRow);
}

// instantiation
template
void cuda_array_math_abs<int32_t >(int32_t*  const array, size_t nRow);
template
void cuda_array_math_abs<uint32_t>(uint32_t* const array, size_t nRow);
template
void cuda_array_math_abs<int64_t >(int64_t*  const array, size_t nRow);
template
void cuda_array_math_abs<uint64_t>(uint64_t* const array, size_t nRow);
template
void cuda_array_math_abs<float   >(float*    const array, size_t nRow);
template
void cuda_array_math_abs<double  >(double*   const array, size_t nRow);
template
void cuda_array_math_abs<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_abs<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_math_abs<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_abs<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_math_abs<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_math_abs<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// inverse
template <typename T> __global__
void dev_array_math_inv(T* const array, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] = -array[iRow];
    }
}

template <typename T> __global__
void dev_array_math_inv(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = -arraySrc[iRow];
    }
}

template <typename T>
void cuda_array_math_inv(T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    dev_array_math_inv
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array, nRow);
}

template <typename T>
void cuda_array_math_inv(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_inv
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, nRow);
}

// instantiation
template
void cuda_array_math_inv<int32_t >(int32_t*  const array, size_t nRow);
template
void cuda_array_math_inv<uint32_t>(uint32_t* const array, size_t nRow);
template
void cuda_array_math_inv<int64_t >(int64_t*  const array, size_t nRow);
template
void cuda_array_math_inv<uint64_t>(uint64_t* const array, size_t nRow);
template
void cuda_array_math_inv<float   >(float*    const array, size_t nRow);
template
void cuda_array_math_inv<double  >(double*   const array, size_t nRow);
template
void cuda_array_math_inv<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_inv<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_math_inv<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_inv<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_math_inv<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_math_inv<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// norm (mean = 0, std = 1)
template <typename T>
void cuda_array_math_norm(T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    T         mean, std;

    mean = cuda_array_reduce_sum(array, nRow) / (T) nRow;

    cuda_array_math_sub(array, mean, nRow);

    std = (T) std::sqrt(cuda_array_reduce_sum2(array, nRow) / (T) nRow);

    cuda_array_math_div(array, std, nRow);
}

template <typename T>
void cuda_array_math_norm(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    T         mean, std;

    mean = cuda_array_reduce_sum(arraySrc, nRow) / (T) nRow;

    cuda_array_math_sub(arrayDst, arraySrc, mean, nRow);

    std = (T) std::sqrt(cuda_array_reduce_sum2(arrayDst, nRow) / (T) nRow);

    cuda_array_math_div(arrayDst, std, nRow);
}

// instantiation
template
void cuda_array_math_norm<int32_t >(int32_t*  const array, size_t nRow);
template
void cuda_array_math_norm<uint32_t>(uint32_t* const array, size_t nRow);
template
void cuda_array_math_norm<int64_t >(int64_t*  const array, size_t nRow);
template
void cuda_array_math_norm<uint64_t>(uint64_t* const array, size_t nRow);
template
void cuda_array_math_norm<float   >(float*    const array, size_t nRow);
template
void cuda_array_math_norm<double  >(double*   const array, size_t nRow);
template
void cuda_array_math_norm<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_norm<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_math_norm<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_norm<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_math_norm<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_math_norm<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// square
template <typename T> __global__
void dev_array_math_sqr(T* const array, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] *= array[iRow];
    }
}

template <typename T> __global__
void dev_array_math_sqr(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = arraySrc[iRow]*arraySrc[iRow];
    }
}

template <typename T>
void cuda_array_math_sqr(T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    dev_array_math_sqr
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array, nRow);
}

template <typename T>
void cuda_array_math_sqr(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_sqr
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, nRow);
}

// instantiation
template
void cuda_array_math_sqr<int32_t >(int32_t*  const array, size_t nRow);
template
void cuda_array_math_sqr<uint32_t>(uint32_t* const array, size_t nRow);
template
void cuda_array_math_sqr<int64_t >(int64_t*  const array, size_t nRow);
template
void cuda_array_math_sqr<uint64_t>(uint64_t* const array, size_t nRow);
template
void cuda_array_math_sqr<float   >(float*    const array, size_t nRow);
template
void cuda_array_math_sqr<double  >(double*   const array, size_t nRow);
template
void cuda_array_math_sqr<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_sqr<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_math_sqr<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_sqr<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_math_sqr<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_math_sqr<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// square root
template <typename T> __global__
void dev_array_math_sqrt(T* const array, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] = sqrtf(array[iRow]);
    }
}

template <> __global__
void dev_array_math_sqrt<double>(double* const array, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] = sqrt(array[iRow]);
    }
}

template <typename T> __global__
void dev_array_math_sqrt(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = sqrtf(arraySrc[iRow]);
    }
}

template <> __global__
void dev_array_math_sqrt<double>(double* const arrayDst, const double* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = sqrt(arraySrc[iRow]);
    }
}

template <typename T>
void cuda_array_math_sqrt(T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    dev_array_math_sqrt
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array, nRow);
}

template <typename T>
void cuda_array_math_sqrt(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_sqrt
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, nRow);
}

// instantiation
template
void cuda_array_math_sqrt<int32_t >(int32_t*  const array, size_t nRow);
template
void cuda_array_math_sqrt<uint32_t>(uint32_t* const array, size_t nRow);
template
void cuda_array_math_sqrt<int64_t >(int64_t*  const array, size_t nRow);
template
void cuda_array_math_sqrt<uint64_t>(uint64_t* const array, size_t nRow);
template
void cuda_array_math_sqrt<float   >(float*    const array, size_t nRow);
template
void cuda_array_math_sqrt<double  >(double*   const array, size_t nRow);
template
void cuda_array_math_sqrt<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_sqrt<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_math_sqrt<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_sqrt<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_math_sqrt<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_math_sqrt<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// add
template <typename T> __global__
void dev_array_math_add(T* const array, T value, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] += value;
    }
}

template <typename T> __global__
void dev_array_math_add(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] += arraySrc[iRow];
    }
}

template <typename T> __global__
void dev_array_math_add(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = arraySrc[iRow] + value;
    }
}

template <typename T> __global__
void dev_array_math_add(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = arraySrc1[iRow] + arraySrc2[iRow];
    }
}

template <typename T>
void cuda_array_math_add(T* const array, T value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    dev_array_math_add
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array, value, nRow);
}

template <typename T>
void cuda_array_math_add(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_add
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, nRow);
}

template <typename T>
void cuda_array_math_add(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_add
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, value, nRow);
}

template <typename T>
void cuda_array_math_add(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(arrayDst  != NULL);
    assert(nRow > 0);

    dev_array_math_add
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc1, arraySrc2, nRow);
}

// instantiation
template
void cuda_array_math_add<int32_t >(int32_t*  const array, int32_t  value, size_t nRow);
template
void cuda_array_math_add<uint32_t>(uint32_t* const array, uint32_t value, size_t nRow);
template
void cuda_array_math_add<int64_t >(int64_t*  const array, int64_t  value, size_t nRow);
template
void cuda_array_math_add<uint64_t>(uint64_t* const array, uint64_t value, size_t nRow);
template
void cuda_array_math_add<float   >(float*    const array, float    value, size_t nRow);
template
void cuda_array_math_add<double  >(double*   const array, double   value, size_t nRow);
template
void cuda_array_math_add<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_add<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_math_add<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_add<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_math_add<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_math_add<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);
template
void cuda_array_math_add<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, int32_t  value, size_t nRow);
template
void cuda_array_math_add<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, uint32_t value, size_t nRow);
template
void cuda_array_math_add<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, int64_t  value, size_t nRow);
template
void cuda_array_math_add<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, uint64_t value, size_t nRow);
template
void cuda_array_math_add<float   >(float*    const arrayDst, const float*    const arraySrc, float    value, size_t nRow);
template
void cuda_array_math_add<double  >(double*   const arrayDst, const double*   const arraySrc, double   value, size_t nRow);
template
void cuda_array_math_add<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow);
template
void cuda_array_math_add<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow);
template
void cuda_array_math_add<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow);
template
void cuda_array_math_add<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow);
template
void cuda_array_math_add<float   >(float*    const arrayDst, const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow);
template
void cuda_array_math_add<double  >(double*   const arrayDst, const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow);

// substract
template <typename T> __global__
void dev_array_math_sub(T* const array, T value, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] -= value;
    }
}

template <typename T> __global__
void dev_array_math_sub(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] -= arraySrc[iRow];
    }
}

template <typename T> __global__
void dev_array_math_sub(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = arraySrc[iRow] - value;
    }
}

template <typename T> __global__
void dev_array_math_sub(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = arraySrc1[iRow] - arraySrc2[iRow];
    }
}

template <typename T>
void cuda_array_math_sub(T* const array, T value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    dev_array_math_sub
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array, value, nRow);
}

template <typename T>
void cuda_array_math_sub(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_sub
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, nRow);
}

template <typename T>
void cuda_array_math_sub(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_sub
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, value, nRow);
}

template <typename T>
void cuda_array_math_sub(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(arrayDst  != NULL);
    assert(nRow > 0);

    dev_array_math_sub
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc1, arraySrc2, nRow);
}

// instantiation
template
void cuda_array_math_sub<int32_t >(int32_t*  const array, int32_t  value, size_t nRow);
template
void cuda_array_math_sub<uint32_t>(uint32_t* const array, uint32_t value, size_t nRow);
template
void cuda_array_math_sub<int64_t >(int64_t*  const array, int64_t  value, size_t nRow);
template
void cuda_array_math_sub<uint64_t>(uint64_t* const array, uint64_t value, size_t nRow);
template
void cuda_array_math_sub<float   >(float*    const array, float    value, size_t nRow);
template
void cuda_array_math_sub<double  >(double*   const array, double   value, size_t nRow);
template
void cuda_array_math_sub<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_sub<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_math_sub<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_sub<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_math_sub<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_math_sub<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);
template
void cuda_array_math_sub<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, int32_t  value, size_t nRow);
template
void cuda_array_math_sub<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, uint32_t value, size_t nRow);
template
void cuda_array_math_sub<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, int64_t  value, size_t nRow);
template
void cuda_array_math_sub<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, uint64_t value, size_t nRow);
template
void cuda_array_math_sub<float   >(float*    const arrayDst, const float*    const arraySrc, float    value, size_t nRow);
template
void cuda_array_math_sub<double  >(double*   const arrayDst, const double*   const arraySrc, double   value, size_t nRow);
template
void cuda_array_math_sub<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow);
template
void cuda_array_math_sub<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow);
template
void cuda_array_math_sub<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow);
template
void cuda_array_math_sub<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow);
template
void cuda_array_math_sub<float   >(float*    const arrayDst, const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow);
template
void cuda_array_math_sub<double  >(double*   const arrayDst, const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow);

// multiply
template <typename T> __global__
void dev_array_math_mul(T* const array, T value, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] *= value;
    }
}

template <typename T> __global__
void dev_array_math_mul(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] *= arraySrc[iRow];
    }
}

template <typename T> __global__
void dev_array_math_mul(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = arraySrc[iRow] * value;
    }
}

template <typename T> __global__
void dev_array_math_mul(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = arraySrc1[iRow] * arraySrc2[iRow];
    }
}

template <typename T>
void cuda_array_math_mul(T* const array, T value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    dev_array_math_mul
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array, value, nRow);
}

template <typename T>
void cuda_array_math_mul(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_mul
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, nRow);
}

template <typename T>
void cuda_array_math_mul(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_mul
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, value, nRow);
}

template <typename T>
void cuda_array_math_mul(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(arrayDst  != NULL);
    assert(nRow > 0);

    dev_array_math_mul
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc1, arraySrc2, nRow);
}

// instantiation
template
void cuda_array_math_mul<int32_t >(int32_t*  const array, int32_t  value, size_t nRow);
template
void cuda_array_math_mul<uint32_t>(uint32_t* const array, uint32_t value, size_t nRow);
template
void cuda_array_math_mul<int64_t >(int64_t*  const array, int64_t  value, size_t nRow);
template
void cuda_array_math_mul<uint64_t>(uint64_t* const array, uint64_t value, size_t nRow);
template
void cuda_array_math_mul<float   >(float*    const array, float    value, size_t nRow);
template
void cuda_array_math_mul<double  >(double*   const array, double   value, size_t nRow);
template
void cuda_array_math_mul<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_mul<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_math_mul<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_mul<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_math_mul<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_math_mul<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);
template
void cuda_array_math_mul<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, int32_t  value, size_t nRow);
template
void cuda_array_math_mul<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, uint32_t value, size_t nRow);
template
void cuda_array_math_mul<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, int64_t  value, size_t nRow);
template
void cuda_array_math_mul<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, uint64_t value, size_t nRow);
template
void cuda_array_math_mul<float   >(float*    const arrayDst, const float*    const arraySrc, float    value, size_t nRow);
template
void cuda_array_math_mul<double  >(double*   const arrayDst, const double*   const arraySrc, double   value, size_t nRow);
template
void cuda_array_math_mul<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow);
template
void cuda_array_math_mul<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow);
template
void cuda_array_math_mul<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow);
template
void cuda_array_math_mul<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow);
template
void cuda_array_math_mul<float   >(float*    const arrayDst, const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow);
template
void cuda_array_math_mul<double  >(double*   const arrayDst, const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow);

// divide
template <typename T> __global__
void dev_array_math_div(T* const array, T value, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] /= value;
    }
}

template <typename T> __global__
void dev_array_math_div(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] /= arraySrc[iRow];
    }
}

template <typename T> __global__
void dev_array_math_div(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = arraySrc[iRow] / value;
    }
}

template <typename T> __global__
void dev_array_math_div(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayDst[iRow] = arraySrc1[iRow] / arraySrc2[iRow];
    }
}

template <typename T>
void cuda_array_math_div(T* const array, T value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    dev_array_math_div
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array, value, nRow);
}

template <typename T>
void cuda_array_math_div(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_div
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, nRow);
}

template <typename T>
void cuda_array_math_div(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_math_div
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, value, nRow);
}

template <typename T>
void cuda_array_math_div(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(arrayDst  != NULL);
    assert(nRow > 0);

    dev_array_math_div
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc1, arraySrc2, nRow);
}

// instantiation
template
void cuda_array_math_div<int32_t >(int32_t*  const array, int32_t  value, size_t nRow);
template
void cuda_array_math_div<uint32_t>(uint32_t* const array, uint32_t value, size_t nRow);
template
void cuda_array_math_div<int64_t >(int64_t*  const array, int64_t  value, size_t nRow);
template
void cuda_array_math_div<uint64_t>(uint64_t* const array, uint64_t value, size_t nRow);
template
void cuda_array_math_div<float   >(float*    const array, float    value, size_t nRow);
template
void cuda_array_math_div<double  >(double*   const array, double   value, size_t nRow);
template
void cuda_array_math_div<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_div<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void cuda_array_math_div<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void cuda_array_math_div<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void cuda_array_math_div<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void cuda_array_math_div<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);
template
void cuda_array_math_div<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, int32_t  value, size_t nRow);
template
void cuda_array_math_div<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, uint32_t value, size_t nRow);
template
void cuda_array_math_div<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, int64_t  value, size_t nRow);
template
void cuda_array_math_div<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, uint64_t value, size_t nRow);
template
void cuda_array_math_div<float   >(float*    const arrayDst, const float*    const arraySrc, float    value, size_t nRow);
template
void cuda_array_math_div<double  >(double*   const arrayDst, const double*   const arraySrc, double   value, size_t nRow);
template
void cuda_array_math_div<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow);
template
void cuda_array_math_div<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow);
template
void cuda_array_math_div<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow);
template
void cuda_array_math_div<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow);
template
void cuda_array_math_div<float   >(float*    const arrayDst, const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow);
template
void cuda_array_math_div<double  >(double*   const arrayDst, const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow);

} // namespace gem

/***********************************************************************
 *  File:       array_random.cu
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

namespace gem {

/*****************************************
 * Random
 ****************************************/

// CURAND setup
template <typename T> __global__
void dev_array_random_setup(T* const globalState, size_t seed)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    curand_init(seed, iRow, 0, &globalState[iRow]);
}

// normal distribution
template <typename T> __global__
void dev_array_random_normal(curandState* const globalState, T* const array, size_t nRow, T mean, T std)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        curandState localState = globalState[iRow];

        array[iRow] = curand_normal(&localState) * std + mean;

        globalState[iRow] = localState;
    }
}

template <> __global__
void dev_array_random_normal<double>(curandState* const globalState, double* const array, size_t nRow, double mean, double std)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        curandState localState = globalState[iRow];

        array[iRow] = curand_normal_double(&localState) * std + mean;

        globalState[iRow] = localState;
    }
}

template <typename T>
void cuda_array_random_normal(T* const array, size_t nRow, T mean, T std, size_t seed)
{
    assert(array != NULL);
    assert(nRow > 0);

    curandState *globalState;

    cuda_arrayDev_new(globalState, nRow);

    dev_array_random_setup
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (globalState, seed);

    dev_array_random_normal
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (globalState, array, nRow, mean, std);

    cuda_arrayDev_delete(globalState);
}

// instantiation
template
void cuda_array_random_normal<float >(float*  const array, size_t nRow, float  mean, float  std, size_t seed);
template
void cuda_array_random_normal<double>(double* const array, size_t nRow, double mean, double std, size_t seed);

// uniform distribution
template <typename T> __global__
void dev_array_random_uniform(curandState* const globalState, T* const array, size_t nRow, T min, T range)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        curandState localState = globalState[iRow];

        array[iRow] = curand_uniform(&localState) * range + min;

        globalState[iRow] = localState;
    }
}

template <> __global__
void dev_array_random_uniform<double>(curandState* const globalState, double* const array, size_t nRow, double min, double range)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        curandState localState = globalState[iRow];

        array[iRow] = curand_uniform_double(&localState) * range + min;

        globalState[iRow] = localState;
    }
}

template <typename T>
void cuda_array_random_uniform(T* const array, size_t nRow, T min, T max, size_t seed)
{
    assert(array != NULL);
    assert(nRow > 0);
    assert(min <= max);

    curandState* globalState = NULL;

    cuda_arrayDev_new(globalState, nRow);

    dev_array_random_setup
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (globalState, seed);

    dev_array_random_uniform
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (globalState, array, nRow, min, max-min);

    cuda_arrayDev_delete(globalState);
}

// instantiation
template
void cuda_array_random_uniform<float >(float*  const array, size_t nRow, float  min, float  max, size_t seed);
template
void cuda_array_random_uniform<double>(double* const array, size_t nRow, double min, double max, size_t seed);

} // namespace gem

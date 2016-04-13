/***********************************************************************
 *  File:     array_reduction_thrust.cu
 *
 *  Purpose:  Implementation of array-related functions
 *
 *  Author:   Thai V. Hoang
 *
 *  Contact:  hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "array.cuh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#pragma GCC diagnostic pop

namespace gem {

/*****************************************
 * Reduction
 ****************************************/

template <typename T>
struct functor_square
{
    __host__ __device__
    T operator()(const T& x) const { return x * x; }
};

template <typename T>
T cuda_array_reduce_thrust(const T* const array, size_t nRow, eReduce reduce)
{
    assert(array != NULL);
    assert(nRow > 0);

    T                       val = 0;
    thrust::device_ptr<T>   array_thrust(const_cast<T*>(array));

    switch (reduce)
    {
        case REDUCE_SUM:
            val = thrust::reduce(array_thrust, array_thrust + nRow,
                                 (T) 0, thrust::plus<T>());
            break;
        case REDUCE_SUM2:
            val = thrust::transform_reduce(array_thrust, array_thrust + nRow,
                                           functor_square<T>(),
                                           (T) 0, thrust::plus<T>());
            break;
        case REDUCE_MIN:
            val = *thrust::min_element(array_thrust, array_thrust + nRow);
            break;
        case REDUCE_MAX:
            val = *thrust::max_element(array_thrust, array_thrust + nRow);
            break;
        case REDUCE_MIN_INDEX:
            val = (T) (thrust::min_element(array_thrust, array_thrust + nRow) - array_thrust);
            break;
        case REDUCE_MAX_INDEX:
            val = (T) (thrust::max_element(array_thrust, array_thrust + nRow) - array_thrust);
            break;
        default:
            ERROR("cuda_array_reduce_thrust", "unsupported reduction mode");
    }

    return val;
}

// instantiation
template
int32_t  cuda_array_reduce_thrust<int32_t >(const int32_t*  const array, size_t nRow, eReduce reduce);
template
uint32_t cuda_array_reduce_thrust<uint32_t>(const uint32_t* const array, size_t nRow, eReduce reduce);
template
int64_t  cuda_array_reduce_thrust<int64_t >(const int64_t*  const array, size_t nRow, eReduce reduce);
template
uint64_t cuda_array_reduce_thrust<uint64_t>(const uint64_t* const array, size_t nRow, eReduce reduce);
template
float    cuda_array_reduce_thrust<float   >(const float*    const array, size_t nRow, eReduce reduce);
template
double   cuda_array_reduce_thrust<double  >(const double*   const array, size_t nRow, eReduce reduce);

} // namespace gem

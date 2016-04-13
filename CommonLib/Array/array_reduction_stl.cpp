/***********************************************************************
 *  File:       array_reduction_stl.cpp
 *
 *  Purpose:    Implementation of array-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "array.hpp"

#include <algorithm>
#include <numeric>
#include <parallel/algorithm>
#include <parallel/numeric>

// Note: Using GNU parallel mode for STL algorithm
//      https://gcc.gnu.org/onlinedocs/libstdc++/manual/parallel_mode.html

namespace gem {

/*****************************************
 * Reduction
 ****************************************/

template <typename T>
struct functor_square
{
    T operator()(const T& left, const T& right) const { return left + right * right; }
};

template <typename T>
T array_reduce_stl(const T* const array, size_t nRow, eReduce reduce, bool parallel)
{
    assert(array != NULL);
    assert(nRow > 0);

    T       val = 0;

    if (parallel) {
        switch (reduce)
        {
            case REDUCE_SUM:
                val = __gnu_parallel::accumulate(array, array + nRow, (T) 0);
                break;
            case REDUCE_SUM2:
                val = __gnu_parallel::accumulate(array, array + nRow, (T) 0, functor_square<T>());
                break;
            case REDUCE_MIN:
                val = *__gnu_parallel::min_element(array, array + nRow);
                break;
            case REDUCE_MAX:
                val = *__gnu_parallel::max_element(array, array + nRow);
                break;
            case REDUCE_MIN_INDEX:
                val = (T) (__gnu_parallel::min_element(array, array + nRow) - array);
                break;
            case REDUCE_MAX_INDEX:
                val = (T) (__gnu_parallel::max_element(array, array + nRow) - array);
                break;
            default:
                ERROR("array_reduce_stl", "unsupported reduction mode");
        }
    }
    else {
        switch (reduce)
        {
            case REDUCE_SUM:
                val = std::accumulate(array, array + nRow, (T) 0);
                break;
            case REDUCE_SUM2:
                val = std::accumulate(array, array + nRow, (T) 0, functor_square<T>());
                break;
            case REDUCE_MIN:
                val = *std::min_element(array, array + nRow);
                break;
            case REDUCE_MAX:
                val = *std::max_element(array, array + nRow);
                break;
            case REDUCE_MIN_INDEX:
                val = (T) (std::min_element(array, array + nRow) - array);
                break;
            case REDUCE_MAX_INDEX:
                val = (T) (std::max_element(array, array + nRow) - array);
                break;
            default:
                ERROR("array_reduce_stl", "unsupported reduction mode");
        }
    }

    return val;
}

// instantiation
template
int32_t  array_reduce_stl<int32_t >(const int32_t*  const array, size_t nRow, eReduce reduce, bool parallel);
template
uint32_t array_reduce_stl<uint32_t>(const uint32_t* const array, size_t nRow, eReduce reduce, bool parallel);
template
int64_t  array_reduce_stl<int64_t >(const int64_t*  const array, size_t nRow, eReduce reduce, bool parallel);
template
uint64_t array_reduce_stl<uint64_t>(const uint64_t* const array, size_t nRow, eReduce reduce, bool parallel);
template
float    array_reduce_stl<float   >(const float*    const array, size_t nRow, eReduce reduce, bool parallel);
template
double   array_reduce_stl<double  >(const double*   const array, size_t nRow, eReduce reduce, bool parallel);

} // namespace gem

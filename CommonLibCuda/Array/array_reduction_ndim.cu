/***********************************************************************
 *  File:       array_reduction_ndim.cu
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
 * Reduction NDim
 ****************************************/

template <typename T>
void cuda_array_reducendim_max (const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim);
template <typename T>
void cuda_array_reducendim_min (const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim);
template <typename T>
void cuda_array_reducendim_sum (const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim);
template <typename T>
void cuda_array_reducendim_sum2(const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim);
template <typename T>
void cuda_array_reducendim_mean(const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim);
template <typename T>
void cuda_array_reducendim_std (const T* const arraySrc, size_t nRow, size_t nCol,
                                      T* const arrayDst,
                                eReduceNDim dim);

template <typename T>
void cuda_array_reducendim_max (const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim);
template <typename T>
void cuda_array_reducendim_min (const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim);
template <typename T>
void cuda_array_reducendim_sum (const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim);
template <typename T>
void cuda_array_reducendim_sum2(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim);
template <typename T>
void cuda_array_reducendim_mean(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim);
template <typename T>
void cuda_array_reducendim_std(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      T* const arrayDst,
                                eReduceNDim dim);

} // namespace gem

/***********************************************************************
 *  File:       array.hpp
 *
 *  Purpose:    Header file for array-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_LIB_ARRAY_HPP__
#define __GEM_LIB_ARRAY_HPP__

#include "config.hpp"
#include "macro.hpp"

#ifdef __GEM_USE_OPENMP__
#include <omp.h>
#endif

#include <stdint.h>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <complex>
#include <string>

namespace gem {

enum ePermute
{
    PERMUTE2D     = 0,
    PERMUTE3D_123 = 1,
    PERMUTE3D_132 = 2,
    PERMUTE3D_213 = 3,
    PERMUTE3D_231 = 4,
    PERMUTE3D_312 = 5,
    PERMUTE3D_321 = 6
};

enum eReduce
{
    REDUCE_SUM       = 0,
    REDUCE_SUM2      = 1,
    REDUCE_MIN       = 2,
    REDUCE_MIN_INDEX = 3,
    REDUCE_MAX       = 4,
    REDUCE_MAX_INDEX = 5,
    REDUCE_MEAN      = 6,
    REDUCE_STD       = 7
};

enum eReduceNDim
{
    REDUCE_NDIM_ROW = 0,
    REDUCE_NDIM_COL = 1,
    REDUCE_NDIM_SEC = 2
};

/*****************************************
 * memcpy / memset / memcheck
 ****************************************/

// memcpy
template <typename T>
void array_memcpy(T* const arrayDst, const T* const arraySrc, size_t nRow);

// memset
template <typename T>
void array_memset(T* const array, int value, size_t nRow);

// check memory leak
void array_memcheck(void);

/*****************************************
 * Array allocation/deallocation
 *       w/o zero initialization
 ****************************************/

// size(array1D) = [nRow]
template <typename T>
void array_new(T*& array1D, size_t nRow);

template <typename T>
void array_new_zero(T*& array1D, size_t nRow);

template <typename T>
void array_delete(T*& array1D);

// size(array2D) = [nRow nCol]
template <typename T>
void array_new(T**& array2D, T*& array1D, size_t nRow, size_t nCol);

template <typename T>
void array_new_zero(T**& array2D, T*& array1D, size_t nRow, size_t nCol);

template <typename T>
void array_delete(T**& array2D, T*& array1D);

// size(array3D) = [nRow nCol nSec]
template <typename T>
void array_new(T***& array3D, T*& array1D, size_t nRow, size_t nCol, size_t nSec);

template <typename T>
void array_new_zero(T***& array3D, T*& array1D, size_t nRow, size_t nCol, size_t nSec);

template <typename T>
void array_delete(T***& array3D, T*& array1D, size_t nRow);

/*****************************************
 * Padding
 ****************************************/

// 1D
template <typename T1, typename T2>
void array_pad(const T1* const arraySrc, size_t nRowSrc,
                     T2* const arrayDst, size_t nRowDst,
                                         size_t nRowOff);

// 2D
template <typename T1, typename T2>
void array_pad(const T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                     T2* const arrayDst, size_t nRowDst, size_t nColDst,
                                         size_t nRowOff, size_t nColOff);

// 3D
template <typename T1, typename T2>
void array_pad(const T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                     T2* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                         size_t nRowOff, size_t nColOff, size_t nSecOff);

/*****************************************
 * Cropping
 ****************************************/

// 1D
template <typename T1, typename T2>
void array_crop(const T1* const arraySrc, size_t nRowSrc,
                      T2* const arrayDst, size_t nRowDst,
                                          size_t nRowOff);

// 2D
template <typename T1, typename T2>
void array_crop(const T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                      T2* const arrayDst, size_t nRowDst, size_t nColDst,
                                          size_t nRowOff, size_t nColOff);

// 3D
template <typename T1, typename T2>
void array_crop(const T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                      T2* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                          size_t nRowOff, size_t nColOff, size_t nSecOff);

/*****************************************
 * Replacement
 ****************************************/

// 1D value
template <typename T1, typename T2>
void array_replace(T1* const arraySrc, size_t nRowSrc,
                   T2        value,    size_t nRowRep,
                                       size_t nRowOff);

// 2D value
template <typename T1, typename T2>
void array_replace(T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                   T2        value,    size_t nRowRep, size_t nColRep,
                                       size_t nRowOff, size_t nColOff);

// 3D value
template <typename T1, typename T2>
void array_replace(T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                   T2        value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                       size_t nRowOff, size_t nColOff, size_t nSecOff);

// 1D array
template <typename T1, typename T2>
void array_replace(      T1* const arraySrc, size_t nRowSrc,
                   const T2* const arrayRep, size_t nRowRep,
                                             size_t nRowOff);

// 2D array
template <typename T1, typename T2>
void array_replace(      T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                   const T2* const arrayRep, size_t nRowRep, size_t nColRep,
                                             size_t nRowOff, size_t nColOff);

// 3D array
template <typename T1, typename T2>
void array_replace(      T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                   const T2* const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                             size_t nRowOff, size_t nColOff, size_t nSecOff);

/*****************************************
 * Indexing
 ****************************************/

// 1D index -> 2D sub-indices
void array_index_ind2sub(const size_t* const arrayInd,
                               size_t* const arrayRow,
                               size_t* const arrayCol,
                               size_t        length,
                               size_t        nRow, size_t nCol);

// 1D index -> 3D sub-indices
void array_index_ind2sub(const size_t* const arrayInd,
                               size_t* const arrayRow,
                               size_t* const arrayCol,
                               size_t* const arraySec,
                               size_t        length,
                               size_t        nRow, size_t nCol, size_t nSec);

// 2D sub-indices -> 1D index
void array_index_sub2ind(const size_t* const arrayRow,
                         const size_t* const arrayCol,
                               size_t* const arrayInd,
                               size_t        length,
                               size_t        nRow, size_t nCol);

// 3D sub-indices -> 1D index
void array_index_sub2ind(const size_t* const arrayRow,
                         const size_t* const arrayCol,
                         const size_t* const arraySec,
                               size_t* const arrayInd,
                               size_t        length,
                               size_t        nRow, size_t nCol, size_t nSec);

// column-major -> row-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow)
template <typename T1, typename T2>
void array_index_col2row(const T1* const arrayCol, T2* const arrayRow,
                         size_t nRow, size_t nCol);

// column-major -> row-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow), nSec = dim3(arrayRow)
template <typename T1, typename T2>
void array_index_col2row(const T1* const arrayCol, T2* const arrayRow,
                         size_t nRow, size_t nCol, size_t nSec);

// row-major -> column-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow)
template <typename T1, typename T2>
void array_index_row2col(const T1* const arrayRow, T2* const arrayCol,
                         size_t nRow, size_t nCol);

// row-major -> column-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow), nSec = dim3(arrayRow)
template <typename T1, typename T2>
void array_index_row2col(const T1* const arrayRow, T2* const arrayCol,
                         size_t nRow, size_t nCol, size_t nSec);

/*****************************************
 * Permutation
 ****************************************/

// 2D
template <typename T>
void array_permute(const T* const arraySrc, T* const arrayDst,
                   size_t nRow, size_t nCol);

// 3D
template <typename T>
void array_permute(const T* const arraySrc, T* const arrayDst,
                   size_t nRow, size_t nCol, size_t nSec,
                   ePermute permute);

/*****************************************
 * Element-wise operations
 ****************************************/

// find
template <typename T>
void array_find(const T*      const arraySrc, size_t nRow,
                      size_t* const arrayDst,
                      size_t&       nFound);

// ismember
template <typename T>
void array_ismember(const T* const arraySrc, size_t nRowSrc,
                    const T* const arraySet, size_t nRowSet,
                          T* const arrayMem);

// set value
template <typename T>
void array_setval(T* const array, T value, size_t nRow);

// type cast
template <typename T1, typename T2>
void array_typecast(const T1* const arraySrc, size_t nRow,
                          T2* const arrayDst);

// threshold elements to 0/1
template <typename T>
void array_threshold(T* const array, size_t nRow,
                     T        cutoff);
template <typename T1, typename T2>
void array_threshold(const T1* const arraySrc, size_t nRow,
                           T1        cutoff,
                           T2* const arrayDst);
template <typename T1, typename T2>
void array_threshold(const T1* const arraySrc, size_t nRow,
                     const T1* const arrayThr,
                           T2* const arrayDst);

// cutoff small elements to value
template <typename T>
void array_cutoff(T* const array, size_t nRow,
                  T        cutoff,
                  T        value);
template <typename T>
void array_cutoff(const T* const arraySrc, size_t nRow,
                        T        cutoff,
                        T* const arrayDst,
                        T        value);
template <typename T>
void array_cutoff(const T* const arraySrc, size_t nRow,
                  const T* const arrayThr,
                        T* const arrayDst,
                        T        value);

// suppress big elements to value
template <typename T>
void array_suppress(T* const array, size_t nRow,
                    T        cutoff,
                    T        value);
template <typename T>
void array_suppress(const T* const arraySrc, size_t nRow,
                          T        cutoff,
                          T* const arrayDst,
                          T        value);
template <typename T>
void array_suppress(const T* const arraySrc, size_t nRow,
                    const T* const arrayThr,
                          T* const arrayDst,
                          T        value);

// remove zero elements by adding noise
template <typename T>
void array_remove_rezo(const T* const arrayNoise, size_t nRow,
                             T* const arrayData,
                             T        cutoff);

/*****************************************
 * Math
 ****************************************/

// abs
template <typename T>
void array_math_abs(T* const array, size_t nRow);
template <typename T>
void array_math_abs(T* const arrayDst, const T* const arraySrc, size_t nRow);

// inverse
template <typename T>
void array_math_inv(T* const array, size_t nRow);
template <typename T>
void array_math_inv(T* const arrayDst, const T* const arraySrc, size_t nRow);

// norm (mean = 0, std = 1)
template <typename T>
void array_math_norm(T* const array, size_t nRow);
template <typename T>
void array_math_norm(T* const arrayDst, const T* const arraySrc, size_t nRow);

// square
template <typename T>
void array_math_sqr(T* const array, size_t nRow);
template <typename T>
void array_math_sqr(T* const arrayDst, const T* const arraySrc, size_t nRow);

// square root
template <typename T>
void array_math_sqrt(T* const array, size_t nRow);
template <typename T>
void array_math_sqrt(T* const arrayDst, const T* const arraySrc, size_t nRow);

// add
template <typename T>
void array_math_add(T* const array, T value, size_t nRow);
template <typename T>
void array_math_add(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void array_math_add(T* const arrayDst, const T* const arraySrc, T value, size_t nRow);
template <typename T>
void array_math_add(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow);

// subtract
template <typename T>
void array_math_sub(T* const array, T value, size_t nRow);
template <typename T>
void array_math_sub(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void array_math_sub(T* const arrayDst, const T* const arraySrc, T value, size_t nRow);
template <typename T>
void array_math_sub(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow);

// multiply
template <typename T>
void array_math_mul(T* const array, T value, size_t nRow);
template <typename T>
void array_math_mul(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void array_math_mul(T* const arrayDst, const T* const arraySrc, T value, size_t nRow);
template <typename T>
void array_math_mul(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow);

// divide
template <typename T>
void array_math_div(T* const array, T value, size_t nRow);
template <typename T>
void array_math_div(T* const arrayDst, const T* const arraySrc, size_t nRow);
template <typename T>
void array_math_div(T* const arrayDst, const T* const arraySrc, T value, size_t nRow);
template <typename T>
void array_math_div(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow);

/*****************************************
 * Linear algebra
 ****************************************/

// matrix multiplication
template <typename T>
void array_linalg_matmul(const T* const arraySrc1, size_t nRowSrc1, size_t nColSrc1,
                         const T* const arraySrc2, size_t nRowSrc2, size_t nColSrc2,
                               T* const arrayDst);

// vectors multiplication to create 2D/3D matrix
template <typename T>
void array_linalg_prodvec2mat(T* const arrayDst, const T* const arrayRow, size_t nRow,
                                                 const T* const arrayCol, size_t nCol);
template <typename T>
void array_linalg_prodvec2mat(T* const arrayDst, const T* const arrayRow, size_t nRow,
                                                 const T* const arrayCol, size_t nCol,
                                                 const T* const arraySec, size_t nSec);

// blas axpy
template <typename T>
void array_linalg_axpy(T* const arrayDst, const T* const arraySrc, T value, size_t nRow);

/*****************************************
 * Reduction
 ****************************************/

// max, min, sum, sum of squares
template <typename T>
T array_reduce_max (const T* const array, size_t nRow);
template <typename T>
T array_reduce_min (const T* const array, size_t nRow);
template <typename T>
T array_reduce_sum (const T* const array, size_t nRow);
template <typename T>
T array_reduce_sum2(const T* const array, size_t nRow);

// mean, std
template <typename T>
T array_reduce_mean(const T* const array, size_t nRow);
template <typename T>
T array_reduce_std (const T* const array, size_t nRow);

// index
template <typename T>
size_t array_reduce_maxindx(const T* const array, size_t nRow);
template <typename T>
size_t array_reduce_minindx(const T* const array, size_t nRow);

// corr, conv
template <typename T>
T array_reduce_corr(const T* const arraySrc1,
                    const T* const arraySrc2,
                    size_t nRow);
template <typename T>
T array_reduce_conv(const T* const arraySrc1,
                    const T* const arraySrc2,
                    size_t nRow);

// compare
template <typename T>
bool array_reduce_compare(const T* const arraySrc1, const T* const arraySrc2, size_t nRow, T threshold);

// stl
template <typename T>
T array_reduce_stl(const T* const array, size_t nRow, eReduce reduce, bool parallel = false);

/*****************************************
 * Reduction NDim
 ****************************************/

template <typename T>
void array_reducendim_max (const T* const arraySrc, size_t nRow, size_t nCol,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_ROW);
template <typename T>
void array_reducendim_min (const T* const arraySrc, size_t nRow, size_t nCol,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_ROW);
template <typename T>
void array_reducendim_sum (const T* const arraySrc, size_t nRow, size_t nCol,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_ROW);
template <typename T>
void array_reducendim_sum2(const T* const arraySrc, size_t nRow, size_t nCol,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_ROW);
template <typename T>
void array_reducendim_mean(const T* const arraySrc, size_t nRow, size_t nCol,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_ROW);
template <typename T>
void array_reducendim_std (const T* const arraySrc, size_t nRow, size_t nCol,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_ROW);

template <typename T>
void array_reducendim_max (const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_SEC);
template <typename T>
void array_reducendim_min (const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_SEC);
template <typename T>
void array_reducendim_sum (const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_SEC);
template <typename T>
void array_reducendim_sum2(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_SEC);
template <typename T>
void array_reducendim_mean(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_SEC);
template <typename T>
void array_reducendim_std(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                 T* const arrayDst,
                           eReduceNDim dim = REDUCE_NDIM_SEC);

/*****************************************
 * Circular shifting
 ****************************************/

// 1D
template <typename T>
void array_circshift(const T* const arraySrc, T* const arrayDst,
                     size_t    nRow,
                     ptrdiff_t nRowOff);

// 2D
template <typename T>
void array_circshift(const T* const arraySrc, T* const arrayDst,
                     size_t    nRow,    size_t    nCol,
                     ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D
template <typename T>
void array_circshift(const T* const arraySrc, T* const arrayDst,
                     size_t    nRow,    size_t    nCol,    size_t    nSec,
                     ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

/*****************************************
 * Shifting
 ****************************************/

// 1D
template <typename T>
void array_shift(const T* const arraySrc, T* const arrayDst,
                 size_t    nRow,
                 ptrdiff_t nRowOff);

// 2D
template <typename T>
void array_shift(const T* const arraySrc, T* const arrayDst,
                 size_t    nRow,    size_t    nCol,
                 ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D
template <typename T>
void array_shift(const T* const arraySrc, T* const arrayDst,
                 size_t    nRow,    size_t    nCol,    size_t    nSec,
                 ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

/*****************************************
 * Mask
 ****************************************/

// 1D copy
template <typename T>
void array_mask_copy(const T* const array,  size_t    nRow,
                     const T* const mask,   size_t    nRowMsk,
                           T* const output, ptrdiff_t nRowOff);

// 2D copy
template <typename T>
void array_mask_copy(const T* const array,  size_t    nRow,    size_t    nCol,
                     const T* const mask,   size_t    nRowMsk, size_t    nColMsk,
                           T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D copy
template <typename T>
void array_mask_copy(const T* const array,  size_t    nRow,    size_t    nCol,    size_t    nSec,
                     const T* const mask,   size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                           T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

// 1D replace
template <typename T>
void array_mask_replace(      T* const array, size_t    nRow,
                        const T* const mask,  size_t    nRowMsk,
                              T        value, ptrdiff_t nRowOff,
                              bool     maskInside = true);

// 2D replace
template <typename T>
void array_mask_replace(      T* const array, size_t    nRow,    size_t    nCol,
                        const T* const mask,  size_t    nRowMsk, size_t    nColMsk,
                              T        value, ptrdiff_t nRowOff, ptrdiff_t nColOff,
                              bool     maskInside = true);
// 3D replace
template <typename T>
void array_mask_replace(      T* const array, size_t    nRow,    size_t    nCol,    size_t    nSec,
                        const T* const mask,  size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                              T        value, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff,
                              bool     maskInside = true);

// sum
template <typename T>
T array_mask_sum(const T* const array, const T* const mask, size_t nRow);

// substract
template <typename T>
void array_mask_sub(T* const arrayDst, const T* const arraySrc, const T* const mask, T value, size_t nRow);

// substract - square - sum
template <typename T>
T array_mask_subsqrsum(const T* const array, T value, const T* const mask, size_t nRow);

// substract - divide - multiply
template <typename T>
void array_mask_subdivmul(T* const array, T valsub, T valdiv, const T* const mask, size_t nRow);

/*****************************************
 * Printing
 ****************************************/

// 1D complex array
template <typename T>
void array_print(const std::complex<T>* const array,
                 size_t nRow,
                 int width = 8,
                 int precision = 4,
                 bool scientific = 0);

// 1D real array
template <typename T>
void array_print(const T* const array,
                 size_t nRow,
                 int width = 8,
                 int precision = 4,
                 bool scientific = 0);

// 2D complex array (*)
template <typename T>
void array_print(const std::complex<T>* const array,
                 size_t nRow, size_t nCol,
                 int width = 8,
                 int precision = 4,
                 bool scientific = 0);

// 2D complex array (**)
template <typename T>
void array_print(std::complex<T> const * const * const array,
                 size_t nRow, size_t nCol,
                 int width = 8,
                 int precision = 4,
                 bool scientific = 0);

// 2D real array (*)
template <typename T>
void array_print(const T* const array,
                 size_t nRow, size_t nCol,
                 int width = 8,
                 int precision = 4,
                 bool scientific = 0);

// 2D real array (**)
template <typename T>
void array_print(T const * const * const array,
                 size_t nRow, size_t nCol,
                 int width = 8,
                 int precision = 4,
                 bool scientific = 0);

// 3D complex array (*)
template <typename T>
void array_print(const std::complex<T>* const array,
                 size_t nRow, size_t nCol, size_t nSec,
                 int width = 8,
                 int precision = 4,
                 bool scientific = 0);

// 3D complex array (***)
template <typename T>
void array_print(std::complex<T> const * const * const * const array,
                 size_t nRow, size_t nCol, size_t nSec,
                 int width = 8,
                 int precision = 4,
                 bool scientific = 0);

// 3D real array (*)
template <typename T>
void array_print(const T* const array,
                 size_t nRow, size_t nCol, size_t nSec,
                 int width = 8,
                 int precision = 4,
                 bool scientific = 0);

// 3D real array (***)
template <typename T>
void array_print(T const * const * const * const array,
                 size_t nRow, size_t nCol, size_t nSec,
                 int width = 8,
                 int precision = 4,
                 bool scientific = 0);

} // namespace gem

#endif

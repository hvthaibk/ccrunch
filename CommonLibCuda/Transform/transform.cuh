/***********************************************************************
 *  File:       transform.cuh
 *
 *  Purpose:    Header file for transformation functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_DEV_TRANSFORM_CUH__
#define __GEM_DEV_TRANSFORM_CUH__

#include "array.cuh"
#include "transform.hpp"

namespace gem {

/**************************
 * Expansion
 *************************/

// 1D
template <typename T>
void cuda_transform_expand_scale_full (const T* const arraySrc, size_t  nRowSrc,
                                             T*&      arrayDst, size_t& nRowDst,
                                             T        factor = 1);

// 2D
template <typename T>
void cuda_transform_expand_rotate_full(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                             T*&      arrayDst, size_t& nRowDst, size_t& nColDst);

template <typename T>
void cuda_transform_expand_scale_full (const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                             T*&      arrayDst, size_t& nRowDst, size_t& nColDst,
                                             T        factor = 1);

// 3D
template <typename T>
void cuda_transform_expand_rotate_full(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                             T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst);

template <typename T>
void cuda_transform_expand_scale_full (const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                             T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst,
                                             T        factor = 1);

/**************************
 * Translation
 *************************/

// 1D
template <typename T>
void cuda_transform_translate(const T* const arraySrc, T* const arrayDst,
                              size_t nRow,
                              T      nRowOff = 0,
                              eInter inter = INTER_LINEAR);

// 2D
template <typename T>
void cuda_transform_translate(const T* const arraySrc, T* const arrayDst,
                              size_t nRow,        size_t nCol,
                              T      nRowOff = 0, T      nColOff = 0,
                              eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void cuda_transform_translate(const T* const arraySrc, T* const arrayDst,
                              size_t nRow,        size_t nCol,        size_t nSec,
                              T      nRowOff = 0, T      nColOff = 0, T      nSecOff = 0,
                              eInter inter = INTER_LINEAR);

/**************************
 * Rotation
 *************************/

// 2D
template <typename T>
void cuda_transform_rotate(const T* const arraySrc, T* const arrayDst,
                           size_t nRow, size_t nCol,
                           T theta = 0,
                           eInter inter = INTER_LINEAR);

void cuda_transform_rotate(const cudaArray* const arraySrc, float* const arrayDst,
                           size_t nRow, size_t nCol,
                           float theta = 0,
                           eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void cuda_transform_rotate(const T* const arraySrc, T* const arrayDst,
                           size_t nRow, size_t nCol, size_t nSec,
                           T alpha = 0, T beta = 0, T gamma = 0,
                           eInter inter = INTER_LINEAR);

void cuda_transform_rotate(const cudaArray* const arraySrc, float* const arrayDst,
                           size_t nRow, size_t nCol, size_t nSec,
                           float alpha = 0, float beta = 0, float gamma = 0,
                           eInter inter = INTER_LINEAR);

/**************************
 * Scaling
 *************************/

// 1D
template <typename T>
void cuda_transform_scale(const T* const arraySrc, T* const arrayDst,
                          size_t nRow,
                          T factor = 1,
                          eInter inter = INTER_LINEAR);

// 2D
template <typename T>
void cuda_transform_scale(const T* const arraySrc, T* const arrayDst,
                          size_t nRow, size_t nCol,
                          T factor = 1,
                          eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void cuda_transform_scale(const T* const arraySrc, T* const arrayDst,
                          size_t nRow, size_t nCol, size_t nSec,
                          T factor = 1,
                          eInter inter = INTER_LINEAR);

/**************************
 * Rotation-Translation
 *************************/

// 2D
template <typename T>
void cuda_transform_rotate_translate(const T* const arraySrc, T* const arrayDst,
                                     size_t nRow,        size_t nCol,
                                     T      theta = 0,
                                     T      nRowOff = 0, T      nColOff = 0,
                                     eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void cuda_transform_rotate_translate(const T* const arraySrc, T* const arrayDst,
                                     size_t nRow,        size_t nCol,        size_t nSec,
                                     T      alpha = 0,   T      beta = 0,    T      gamma = 0,
                                     T      nRowOff = 0, T      nColOff = 0, T      nSecOff = 0,
                                     eInter inter = INTER_LINEAR);

/**************************
 * Interp
 *************************/

// 1D
template <typename T>
void cuda_transform_interp(const T* const arrayRow,
                           const T* const arraySrc,
                                 T* const arrayDst,
                           size_t nRow,
                           eInter inter = INTER_LINEAR);

// 2D
template <typename T>
void cuda_transform_interp(const T* const arrayRow,
                           const T* const arrayCol,
                           const T* const arraySrc,
                                 T* const arrayDst,
                           size_t nRow, size_t nCol,
                           eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void cuda_transform_interp(const T* const arrayRow,
                           const T* const arrayCol,
                           const T* const arraySec,
                           const T* const arraySrc,
                                 T* const arrayDst,
                           size_t nRow, size_t nCol, size_t nSec,
                           eInter inter = INTER_LINEAR);

} // namespace gem

#endif

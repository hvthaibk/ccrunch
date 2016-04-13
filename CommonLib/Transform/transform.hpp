/***********************************************************************
 *  File:       transform.hpp
 *
 *  Purpose:    Header file for transformation functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_LIB_TRANSFORM_HPP__
#define __GEM_LIB_TRANSFORM_HPP__

#include "array.hpp"
#include <set>
#include <vector>

namespace gem {

enum eExpRot
{
    EXP_ROT_VALID = 0,
    EXP_ROT_FULL  = 1
};

enum eInter
{
    INTER_NEAREST = 0,
    INTER_LINEAR  = 1,
    INTER_CUBIC   = 2
};

enum eRot2D
{
    ROT2D_RIGHT = 0
};

enum eRot3D
{
    ROT3D_RIGHT_ZYZ = 0
};

enum eWaterDir
{
    WATER_DIR_RISE = 0,
    WATER_DIR_FALL = 1
};

enum eTrace
{
    TRACE_MAX  = 0,
    TRACE_MIN  = 1,
    TRACE_SUM  = 2,
    TRACE_SUM2 = 3,
    TRACE_MEAN = 4,
    TRACE_STD  = 5
};

/**************************
 * Helper functions
 *************************/

// 2D rotation matrix
template <typename T>
void transform_rotmat(T theta,
                      T& m11, T& m12,
                      T& m21, T& m22,
                      eRot2D order = ROT2D_RIGHT);

// 3D rotation matrix
template <typename T>
void transform_rotmat(T alpha, T beta, T gamma,
                      T& m11,  T& m12, T& m13,
                      T& m21,  T& m22, T& m23,
                      T& m31,  T& m32, T& m33,
                      eRot3D order = ROT3D_RIGHT_ZYZ);

// center of coordinate
template <typename T>
T transform_centerCoord(size_t nRow);

template <typename T>
T transform_centerCoord(T cSrc, T factor);

// scaled size
template <typename T>
size_t transform_scaleSize(size_t nPixSrc, T factor);

template <typename T>
size_t transform_scaleSize(size_t nPixSrc, T cSrc, T factor);

// binned size
size_t transform_binSize(size_t nPixSrc, size_t factor);

/**************************
 * Expansion
 *************************/

// 1D
template <typename T>
void transform_expand_scale_full  (const T* const arraySrc, size_t  nRowSrc,
                                         T*&      arrayDst, size_t& nRowDst,
                                         T        factor = 1);

// 2D
template <typename T>
void transform_expand_rotate_valid(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                         T*&      arrayDst, size_t& nRowDst, size_t& nColDst);

template <typename T>
void transform_expand_rotate_full (const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                         T*&      arrayDst, size_t& nRowDst, size_t& nColDst);

template <typename T>
void transform_expand_scale_full  (const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                         T*&      arrayDst, size_t& nRowDst, size_t& nColDst,
                                         T        factor = 1);

// 3D
template <typename T>
void transform_expand_rotate_valid(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                         T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst);

template <typename T>
void transform_expand_rotate_full (const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                         T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst);

template <typename T>
void transform_expand_scale_full  (const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                         T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst,
                                         T        factor = 1);

/**************************
 * Translation
 *************************/

// 1D
template <typename T>
void transform_translate(const T* const arraySrc, T* const arrayDst,
                         size_t nRow,
                         T      nRowOff = 0,
                         eInter inter = INTER_LINEAR);

// 2D
template <typename T>
void transform_translate(const T* const arraySrc, T* const arrayDst,
                         size_t nRow,        size_t nCol,
                         T      nRowOff = 0, T      nColOff = 0,
                         eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void transform_translate(const T* const arraySrc, T* const arrayDst,
                         size_t nRow,        size_t nCol,        size_t nSec,
                         T      nRowOff = 0, T      nColOff = 0, T      nSecOff = 0,
                         eInter inter = INTER_LINEAR);

/**************************
 * Rotation
 *************************/

// 2D
template <typename T>
void transform_rotate(const T* const arraySrc, T* const arrayDst,
                      size_t nRow, size_t nCol,
                      T theta = 0,
                      eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void transform_rotate(const T* const arraySrc, T* const arrayDst,
                      size_t nRow, size_t nCol, size_t nSec,
                      T alpha = 0, T beta = 0, T gamma = 0,
                      eInter inter = INTER_LINEAR);

/**************************
 * Scaling
 *************************/

// 1D
template <typename T>
void transform_scale(const T* const arraySrc, T* const arrayDst,
                     size_t nRow,
                     T factor = 1,
                     eInter inter = INTER_LINEAR);

template <typename T>
void transform_scale(const T* const arraySrc, T* const arrayDst,
                     size_t nRow,
                     T      cRow,
                     T factor = 1,
                     eInter inter = INTER_LINEAR);

// 2D
template <typename T>
void transform_scale(const T* const arraySrc, T* const arrayDst,
                     size_t nRow, size_t nCol,
                     T factor = 1,
                     eInter inter = INTER_LINEAR);

template <typename T>
void transform_scale(const T* const arraySrc, T* const arrayDst,
                     size_t nRow, size_t nCol,
                     T      cRow, T      cCol,
                     T factor = 1,
                     eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void transform_scale(const T* const arraySrc, T* const arrayDst,
                     size_t nRow, size_t nCol, size_t nSec,
                     T factor = 1,
                     eInter inter = INTER_LINEAR);

template <typename T>
void transform_scale(const T* const arraySrc, T* const arrayDst,
                     size_t nRow, size_t nCol, size_t nSec,
                     T      cRow, T      cCol, T      cSec,
                     T factor = 1,
                     eInter inter = INTER_LINEAR);

/**************************
 * Rotation-Translation
 *************************/

// 2D
template <typename T>
void transform_rotate_translate(const T* const arraySrc, T* const arrayDst,
                                size_t nRow,        size_t nCol,
                                T      theta = 0,
                                T      nRowOff = 0, T      nColOff = 0,
                                eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void transform_rotate_translate(const T* const arraySrc, T* const arrayDst,
                                size_t nRow,        size_t nCol,        size_t nSec,
                                T      alpha = 0,   T      beta = 0,    T      gamma = 0,
                                T      nRowOff = 0, T      nColOff = 0, T      nSecOff = 0,
                                eInter inter = INTER_LINEAR);

/**************************
 * Binning
 *************************/

// 1D
template <typename T>
void transform_bin(const T* const arraySrc, T* const arrayDst,
                   size_t nRow,
                   size_t factor);

// 2D
template <typename T>
void transform_bin(const T* const arraySrc, T* const arrayDst,
                   size_t nRow, size_t nCol,
                   size_t factor);

// 3D
template <typename T>
void transform_bin(const T* const arraySrc, T* const arrayDst,
                   size_t nRow, size_t nCol, size_t nSec,
                   size_t factor);

/**************************
 * Watershed
 *************************/

// 1D/2D/3D
template <typename T>
void transform_watershed(const T*      const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                               size_t* const arrayLbl,
                               size_t        nDim,
                               T             threshold,
                               T             seedDistMin,
                               eWaterDir     water,
                               std::vector< std::set<size_t> >& neighborList);

/**************************
 * Dilation
 *************************/

// 1D
template <typename T>
void transform_dilate(const T* const arraySrc, size_t nRow,
                            T* const arrayDst,
                            T        radius);

// 2D
template <typename T>
void transform_dilate(const T* const arraySrc, size_t nRow, size_t nCol,
                            T* const arrayDst,
                            T        radius);

// 3D
template <typename T>
void transform_dilate(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                            T* const arrayDst,
                            T        radius);

/**************************
 * Labeling
 *************************/

// 1D
template <typename T>
void transform_label(const T*      const arraySrc, size_t nRow,
                           size_t* const labelDst);

// 2D
template <typename T>
void transform_label(const T*      const arraySrc, size_t nRow, size_t nCol,
                           size_t* const labelDst);

// 3D
template <typename T>
void transform_label(const T*      const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                           size_t* const labelDst);

/**************************
 * Distance / Voronoi
 *************************/

// 1D
template <typename T>
void transform_distance(const T*      const arraySrc, size_t nRow,
                              T*      const arrayDst,
                        const size_t* const labelSrc = NULL,
                              size_t* const labelDst = NULL);

// 2D
template <typename T>
void transform_distance(const T*      const arraySrc, size_t nRow, size_t nCol,
                              T*      const arrayDst,
                        const size_t* const labelSrc = NULL,
                              size_t* const labelDst = NULL);

// 3D
template <typename T>
void transform_distance(const T*      const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                              T*      const arrayDst,
                        const size_t* const labelSrc = NULL,
                              size_t* const labelDst = NULL);

/**************************
 * Projection
 *************************/

// 2D
template <typename T>
void transform_project(const T* const arraySrc, size_t nRowSrc, size_t nColSrc,
                             T* const arrayDst,
                       const std::vector<T> arrayTheta,
                       eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void transform_project(const T* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                             T* const arrayDst,
                       const std::vector<T> arrayAlpha,
                       const std::vector<T> arrayBeta,
                       eInter inter = INTER_LINEAR);

// 2D
template <typename T>
void transform_project(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                             T*&      arrayDst, size_t& nRowDst, size_t& nColDst,
                       const std::vector<T> arrayTheta,
                       eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void transform_project(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                             T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst,
                       const std::vector<T> arrayAlpha,
                       const std::vector<T> arrayBeta,
                       eInter inter = INTER_LINEAR);

/**************************
 * Tracing
 *************************/

// 2D
template <typename T>
void transform_trace(const T* const arraySrc, size_t nRowSrc, size_t nColSrc,
                           T* const arrayDst,
                     const std::vector<T> arrayTheta,
                     eTrace trace = TRACE_SUM,
                     eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void transform_trace(const T* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                           T* const arrayDst,
                     const std::vector<T> arrayAlpha,
                     const std::vector<T> arrayBeta,
                     const std::vector<T> arrayGamma,
                     eTrace trace = TRACE_SUM,
                     eInter inter = INTER_LINEAR);

/**************************
 * Interp
 *************************/

// 1D
template <typename T>
void transform_interp(const T* const arrayRow,
                      const T* const arraySrc,
                            T* const arrayDst,
                      size_t nRow,
                      eInter inter = INTER_LINEAR);

// 2D
template <typename T>
void transform_interp(const T* const arrayRow,
                      const T* const arrayCol,
                      const T* const arraySrc,
                            T* const arrayDst,
                      size_t nRow, size_t nCol,
                      eInter inter = INTER_LINEAR);

// 3D
template <typename T>
void transform_interp(const T* const arrayRow,
                      const T* const arrayCol,
                      const T* const arraySec,
                      const T* const arraySrc,
                            T* const arrayDst,
                      size_t nRow, size_t nCol, size_t nSec,
                      eInter inter = INTER_LINEAR);

} // namespace gem

#endif

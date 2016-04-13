/***********************************************************************
 *  File:       basis.cuh
 *
 *  Purpose:    Header file for basis-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_LIB_BASIS_CUH__
#define __GEM_LIB_BASIS_CUH__

#include "array.cuh"

namespace gem {

/**************************
 * Cartesian - Polar
 *************************/

// cart2pol
template <typename T>
void cuda_coord_cart2pol(const T* const arrayX,   const T* const arrayY,
                         T* const arrayRad, T* const arrayThe,
                         size_t length);

template <typename T>
void cuda_coord_cart2pol(const T* const arrayRow, const T* const arrayCol,
                         T* const arrayRad, T* const arrayThe,
                         size_t length,
                         T      cRow, T      cCol,
                         size_t nRow, size_t nCol);

// pol2cart
template <typename T>
void cuda_coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                         T* const arrayX,   T* const arrayY,
                         size_t length);

template <typename T>
void cuda_coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                         T* const arrayRow, T* const arrayCol,
                         size_t length,
                         T      cRow, T      cCol,
                         size_t nRow, size_t nCol);

/**************************
 * Cartesian - Spherical
 *************************/

// cart2sph
template <typename T>
void cuda_coord_cart2sph(const T* const arrayX,   const T* const arrayY,   const T* const arrayZ,
                               T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                         size_t length);

template <typename T>
void cuda_coord_cart2sph(const T* const arrayRow, const T* const arrayCol, const T* const arraySec,
                               T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                         size_t length,
                         T      cRow, T      cCol, T      cSec,
                         size_t nRow, size_t nCol, size_t nSec);

// sph2cart
template <typename T>
void cuda_coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                               T* const arrayX,         T* const arrayY,         T* const arrayZ,
                         size_t length);

template <typename T>
void cuda_coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                               T* const arrayRow,       T* const arrayCol,       T* const arraySec,
                         size_t length,
                         T      cRow, T      cCol, T      cSec,
                         size_t nRow, size_t nCol, size_t nSec);

} // namespace gem

#endif

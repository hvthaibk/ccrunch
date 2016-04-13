/***********************************************************************
 *  File:       basis.hpp
 *
 *  Purpose:    Header file for basis-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_LIB_BASIS_HPP__
#define __GEM_LIB_BASIS_HPP__

#include "array.hpp"

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace gem {

typedef boost::multiprecision::number< boost::multiprecision::cpp_dec_float<100> > floatmp;

enum eBasisAccuracy
{
    BASIS_ACCU_APPROX    = 0,
    BASIS_ACCU_ACCURATE  = 1,
    BASIS_ACCU_MP        = 2
};

enum eBasisFunction
{
    BASIS_CART_GAUSSIAN       = 0,
    BASIS_CART_HERMITE        = 1,
    BASIS_POLAR_RN            = 2,
    BASIS_POLAR_ZM            = 3,
    BASIS_POLAR_PZM           = 4,
    BASIS_SPHER_GAUSSLAGUERRE = 5
};

/*****************************************
 * Basis conversion
 ****************************************/

template <typename T>
void basis_convert_G2H1D(unsigned int nmax, T sigma,
                         T* const arrayCoef);

void basis_convert_G2H1D(unsigned int nmax, floatmp sigma,
                         floatmp* const arrayCoef);

/*****************************************
 * Cartesian basis
 ****************************************/

// 1D
template <typename T>
void basis_cart_gaussian(const T* const arrayCoord, size_t length,
                               T* const arrayBasis,
                         unsigned int nmax, T sigma,
                         eBasisAccuracy accuracy = BASIS_ACCU_MP);

template <typename T>
void basis_cart_hermite (const T* const arrayCoord, size_t length,
                               T* const arrayBasis,
                         unsigned int nmax, T sigma,
                         eBasisAccuracy accuracy = BASIS_ACCU_MP);

void basis_cart_gaussian(const floatmp* const arrayCoord, size_t length,
                               floatmp* const arrayBasis,
                         unsigned int nmax, floatmp sigma);

void basis_cart_hermite (const floatmp* const arrayCoord, size_t length,
                               floatmp* const arrayBasis,
                         unsigned int nmax, floatmp sigma);

/*****************************************
 * Polar basis
 ****************************************/

template <typename T>
void basis_polar_Rn(const T* const arrayCoord,  size_t       length,
                          T* const arrayRadial, unsigned int nmax,
                          T sigma,
                    eBasisFunction method);

template <typename T>
void basis_spher_Am(const T*         const arrayCoord,   size_t       length,
                    std::complex<T>* const arrayAngular, unsigned int nmax);

/*****************************************
 * Spherical basis
 ****************************************/

template <typename T>
void basis_spher_Rln(const T* const arrayCoord,  size_t       length,
                           T* const arrayRadial, unsigned int nmax,
                           T sigma,
                     eBasisFunction method);

template <typename T>
void basis_spher_Ylm(const T*         const arrayCoord,   size_t       length,
                     std::complex<T>* const arrayAngular, unsigned int nmax);

template <typename T>
void basis_spher_ylm(const T* const arrayCoord,   size_t       length,
                           T* const arrayAngular, unsigned int nmax);

/*****************************************
 * Conversion
 ****************************************/

template <typename T>
void basis_convert_gaussian2hermite(void);

/*****************************************
 * Cartesian - Polar
 ****************************************/

// cart2pol
template <typename T>
void coord_cart2pol(const T* const arrayX,   const T* const arrayY,
                          T* const arrayRad,       T* const arrayThe,
                    size_t length);

template <typename T>
void coord_cart2pol(const T* const arrayRow, const T* const arrayCol,
                          T* const arrayRad,       T* const arrayThe,
                    size_t length,
                    T      cRow, T      cCol,
                    size_t nRow, size_t nCol);

// pol2cart
template <typename T>
void coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                          T* const arrayX,         T* const arrayY,
                    size_t length);

template <typename T>
void coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                          T* const arrayRow,       T* const arrayCol,
                    size_t length,
                    T      cRow, T      cCol,
                    size_t nRow, size_t nCol);

/*****************************************
 * Cartesian - Spherical
 *****************************************/

// cart2sph
template <typename T>
void coord_cart2sph(const T* const arrayX,   const  T* const arrayY,   const T* const arrayZ,
                          T* const arrayRad,        T* const arrayThe,       T* const arrayPhi,
                    size_t length);

template <typename T>
void coord_cart2sph(const T* const arrayRow, const T* const arrayCol, const T* const arraySec,
                          T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                    size_t length,
                    T      cRow, T      cCol, T      cSec,
                    size_t nRow, size_t nCol, size_t nSec);

// sph2cart
template <typename T>
void coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                          T* const arrayX,         T* const arrayY,         T* const arrayZ,
                    size_t length);

template <typename T>
void coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                          T* const arrayRow,       T* const arrayCol,       T* const arraySec,
                    size_t length,
                    T      cRow, T      cCol, T      cSec,
                    size_t nRow, size_t nCol, size_t nSec);

} // namespace gem

#endif

/***********************************************************************
 *  File:       coord_spherical.cpp
 *
 *  Purpose:    Implementation of basis-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "basis.hpp"

namespace gem {

/*****************************************
 * Cartesian - Spherical
 *****************************************/

// Reference:
// http://en.wikipedia.org/wiki/Spherical_coordinate_system
// (use PHYSICS notation for compatibility with spherical harmonics)

// cart2sph
template <typename T>
void coord_cart2sph(const T* const arrayX,   const T* const arrayY,   const T* const arrayZ,
                          T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                    size_t length)
{
    assert(arrayX   != NULL && arrayY   != NULL && arrayZ   != NULL);
    assert(arrayRad != NULL && arrayThe != NULL && arrayPhi != NULL);
    assert(length > 0);

    #pragma omp parallel
    for (size_t i = 0; i < length; i++) {
        arrayRad[i] = std::sqrt(pow2(arrayX[i])+pow2(arrayY[i])+pow2(arrayZ[i]));
        arrayThe[i] = std::atan2(std::sqrt(pow2(arrayX[i])+pow2(arrayY[i])),arrayZ[i]);
        arrayPhi[i] = std::atan2(arrayY[i],arrayX[i]);
    }
}

#ifndef NDEBUG
template <typename T>
void coord_cart2sph(const T* const arrayRow, const T* const arrayCol, const T* const arraySec,
                          T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                    size_t length,
                    T      cRow, T      cCol, T      cSec,
                    size_t nRow, size_t nCol, size_t nSec)
#else
template <typename T>
void coord_cart2sph(const T* const arrayRow, const T* const arrayCol, const T* const arraySec,
                          T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                    size_t length,
                    T      cRow, T      cCol, T      cSec,
                    size_t , size_t , size_t )
#endif
{
    assert(arrayRow != NULL && arrayCol != NULL && arraySec != NULL);
    assert(arrayRad != NULL && arrayThe != NULL && arrayPhi != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    T    x, y, z;

    #pragma omp parallel for private(x,y,z)
    for (size_t i = 0; i < length; i++) {
        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(arraySec[i] >= 0 && arraySec[i] < nSec);
        assert(arrayRow[i] == round(arrayRow[i]));
        assert(arrayCol[i] == round(arrayCol[i]));
        assert(arraySec[i] == round(arraySec[i]));

        x = cSec - arraySec[i];
        y = arrayCol[i] - cCol;
        z = cRow - arrayRow[i];

        arrayRad[i] = std::sqrt(pow2(x)+pow2(y)+pow2(z));
        arrayThe[i] = std::atan2(std::sqrt(pow2(x)+pow2(y)),z);
        arrayPhi[i] = std::atan2(y,x);
    }
}

// instantiation
template
void coord_cart2sph<float >(const float*  const arrayRow, const float*  const arrayCol, const float*  const arraySec,
                                  float*  const arrayRad,       float*  const arrayThe,       float*  const arrayPhi,
                            size_t length);
template
void coord_cart2sph<double>(const double* const arrayRow, const double* const arrayCol, const double* const arraySec,
                                  double* const arrayRad,       double* const arrayThe,       double* const arrayPhi,
                            size_t length);
template
void coord_cart2sph<float >(const float*  const arrayRow, const float*  const arrayCol, const float*  const arraySec,
                                  float*  const arrayRad,       float*  const arrayThe,       float*  const arrayPhi,
                            size_t length,
                            float  cRow, float  cCol, float  cSec,
                            size_t nRow, size_t nCol, size_t nSec);
template
void coord_cart2sph<double>(const double* const arrayRow, const double* const arrayCol, const double* const arraySec,
                                  double* const arrayRad,       double* const arrayThe,       double* const arrayPhi,
                            size_t length,
                            double cRow, double cCol, double cSec,
                            size_t nRow, size_t nCol, size_t nSec);

// sph2cart
template <typename T>
void coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                          T* const arrayX,         T* const arrayY,         T* const arrayZ,
                    size_t length)
{
    assert(arrayRad != NULL && arrayThe != NULL && arrayPhi != NULL);
    assert(arrayX   != NULL && arrayY   != NULL && arrayZ   != NULL);
    assert(length > 0);

    T    proj;

    #pragma omp parallel for private(proj)
    for (size_t i = 0; i < length; i++) {
        assert(arrayRad[i] >= 0);

        proj = arrayRad[i]*std::sin(arrayThe[i]);

        arrayX[i] = proj*std::cos(arrayPhi[i]);
        arrayY[i] = proj*std::sin(arrayPhi[i]);
        arrayZ[i] = arrayRad[i]*std::cos(arrayThe[i]);
    }
}

#ifndef NDEBUG
template <typename T>
void coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                          T* const arrayRow,       T* const arrayCol,       T* const arraySec,
                    size_t length,
                    T      cRow, T      cCol, T      cSec,
                    size_t nRow, size_t nCol, size_t nSec)
#else
template <typename T>
void coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                          T* const arrayRow,       T* const arrayCol,       T* const arraySec,
                    size_t length,
                    T      cRow, T      cCol, T      cSec,
                    size_t , size_t , size_t )
#endif
{
    assert(arrayRad != NULL && arrayThe != NULL && arrayPhi != NULL);
    assert(arrayRow != NULL && arrayCol != NULL && arraySec != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    T    proj, tmpRow, tmpCol, tmpSec;

    #pragma omp parallel for private(proj,tmpRow,tmpCol,tmpSec)
    for (size_t i = 0; i < length; i++) {
        assert(arrayRad[i] >= 0);

        proj = arrayRad[i]*std::sin(arrayThe[i]);

        tmpRow = cRow - arrayRad[i]*std::cos(arrayThe[i]);
        tmpCol = cCol + proj*std::sin(arrayPhi[i]);
        tmpSec = cSec - proj*std::cos(arrayPhi[i]);
        arrayRow[i] = round(tmpRow);
        arrayCol[i] = round(tmpCol);
        arraySec[i] = round(tmpSec);

        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(arraySec[i] >= 0 && arraySec[i] < nSec);
        assert(std::abs(arrayRow[i]-tmpRow) < 0.001);
        assert(std::abs(arrayCol[i]-tmpCol) < 0.001);
        assert(std::abs(arraySec[i]-tmpSec) < 0.001);
    }
}

// instantiation
template
void coord_sph2cart<float >(const float*  const arrayRad, const float*  const arrayThe, const float*  const arrayPhi,
                                  float*  const arrayRow,       float*  const arrayCol,       float*  const arraySec,
                            size_t length);
template
void coord_sph2cart<double>(const double* const arrayRad, const double* const arrayThe, const double* const arrayPhi,
                                  double* const arrayRow,       double* const arrayCol,       double* const arraySec,
                            size_t length);
template
void coord_sph2cart<float >(const float*  const arrayRad, const float*  const arrayThe, const float*  const arrayPhi,
                                  float*  const arrayRow,       float*  const arrayCol,       float*  const arraySec,
                            size_t length,
                            float  cRow, float  cCol, float  cSec,
                            size_t nRow, size_t nCol, size_t nSec);
template
void coord_sph2cart<double>(const double* const arrayRad, const double* const arrayThe, const double* const arrayPhi,
                                  double* const arrayRow,       double* const arrayCol,       double* const arraySec,
                            size_t length,
                            double cRow, double cCol, double cSec,
                            size_t nRow, size_t nCol, size_t nSec);

} // namespace gem

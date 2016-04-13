/***********************************************************************
 *  File:       coord_polar.cpp
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
 * Cartesian - Polar
 ****************************************/

// Reference:
// http://en.wikipedia.org/wiki/Polar_coordinate_system

// cart2pol
template <typename T>
void coord_cart2pol(const T* const arrayX,   const T* const arrayY,
                          T* const arrayRad,       T* const arrayThe,
                    size_t length)
{
    assert(arrayX   != NULL && arrayY   != NULL);
    assert(arrayRad != NULL && arrayThe != NULL);
    assert(length > 0);

    #pragma omp parallel for
    for (size_t i = 0; i < length; i++) {
        arrayRad[i] = std::sqrt(pow2(arrayX[i])+pow2(arrayY[i]));
        arrayThe[i] = std::atan2(arrayY[i],arrayX[i]);
    }
}

#ifndef NDEBUG
template <typename T>
void coord_cart2pol(const T* const arrayRow, const T* const arrayCol,
                          T* const arrayRad,       T* const arrayThe,
                    size_t length,
                    T      cRow, T      cCol,
                    size_t nRow, size_t nCol)
#else
template <typename T>
void coord_cart2pol(const T* const arrayRow, const T* const arrayCol,
                          T* const arrayRad,       T* const arrayThe,
                    size_t length,
                    T      cRow, T      cCol,
                    size_t , size_t )
#endif
{
    assert(arrayRow != NULL && arrayCol != NULL);
    assert(arrayRad != NULL && arrayThe != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0);

    T    x, y;

    #pragma omp parallel for private(x,y)
    for (size_t i = 0; i < length; i++) {
        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(arrayRow[i] == round(arrayRow[i]));
        assert(arrayCol[i] == round(arrayCol[i]));

        x = arrayCol[i] - cCol;
        y = cRow - arrayRow[i];

        arrayRad[i] = std::sqrt(pow2(x)+pow2(y));
        arrayThe[i] = std::atan2(y,x);
    }
}

// instantiation
template
void coord_cart2pol<float >(const float*  const arrayRow, const float*  const arrayCol,
                                  float*  const arrayRad,       float*  const arrayThe,
                            size_t length);
template
void coord_cart2pol<double>(const double* const arrayRow, const double* const arrayCol,
                                  double* const arrayRad,       double* const arrayThe,
                            size_t length);
template
void coord_cart2pol<float >(const float*  const arrayRow, const float*  const arrayCol,
                                  float*  const arrayRad,       float*  const arrayThe,
                            size_t length,
                            float  cRow, float  cCol,
                            size_t nRow, size_t nCol);
template
void coord_cart2pol<double>(const double* const arrayRow, const double* const arrayCol,
                                  double* const arrayRad,       double* const arrayThe,
                            size_t length,
                            double cRow, double cCol,
                            size_t nRow, size_t nCol);

// pol2cart
template <typename T>
void coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                          T* const arrayX,         T* const arrayY,
                    size_t length)
{
    assert(arrayRad != NULL && arrayThe != NULL);
    assert(arrayX   != NULL && arrayY   != NULL);
    assert(length > 0);

    #pragma omp parallel for
    for (size_t i = 0; i < length; i++) {
        assert(arrayRad[i] >= 0);

        arrayX[i] = arrayRad[i]*std::cos(arrayThe[i]);
        arrayY[i] = arrayRad[i]*std::sin(arrayThe[i]);
    }
}

#ifndef NDEBUG
template <typename T>
void coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                          T* const arrayRow,       T* const arrayCol,
                    size_t length,
                    T      cRow, T      cCol,
                    size_t nRow, size_t nCol)
#else
template <typename T>
void coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                          T* const arrayRow,       T* const arrayCol,
                    size_t length,
                    T      cRow, T      cCol,
                    size_t , size_t )
#endif
{
    assert(arrayRad != NULL && arrayThe != NULL);
    assert(arrayRow != NULL && arrayCol != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0);

    T   tmpRow, tmpCol;

    #pragma omp parallel for private(tmpRow,tmpCol)
    for (size_t i = 0; i < length; i++) {
        assert(arrayRad[i] >= 0);

        tmpRow = cRow - arrayRad[i]*std::sin(arrayThe[i]);
        tmpCol = cCol + arrayRad[i]*std::cos(arrayThe[i]);
        arrayRow[i] = round(tmpRow);
        arrayCol[i] = round(tmpCol);

        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(std::abs(arrayRow[i]-tmpRow) < 0.001);
        assert(std::abs(arrayCol[i]-tmpCol) < 0.001);
    }
}

// instantiation
template
void coord_pol2cart<float >(const float*  const arrayRad, const float*  const arrayThe,
                                  float*  const arrayRow,       float*  const arrayCol,
                            size_t length);
template
void coord_pol2cart<double>(const double* const arrayRad, const double* const arrayThe,
                                  double* const arrayRow,       double* const arrayCol,
                            size_t length);
template
void coord_pol2cart<float >(const float*  const arrayRad, const float*  const arrayThe,
                                  float*  const arrayRow,       float*  const arrayCol,
                            size_t length,
                            float  cRow, float  cCol,
                            size_t nRow, size_t nCol);
template
void coord_pol2cart<double>(const double* const arrayRad, const double* const arrayThe,
                                  double* const arrayRow,       double* const arrayCol,
                            size_t length,
                            double cRow, double cCol,
                            size_t nRow, size_t nCol);

} // namespace gem

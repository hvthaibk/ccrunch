/***********************************************************************
 *  File:       basis_polar.cu
 *
 *  Purpose:    Implementation of basis-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "basis.cuh"

namespace gem {

/**************************
 * Cartesian - Polar
 *************************/

// cart2pol
template <typename T> __global__
void dev_coord_cart2pol(const T* const arrayX,   const T* const arrayY,
                              T* const arrayRad,       T* const arrayThe,
                        size_t length)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < length) {
        arrayRad[i] = sqrtf(arrayX[i]*arrayX[i] + arrayY[i]*arrayY[i]);
        arrayThe[i] = atan2f(arrayY[i],arrayX[i]);
    }
}

template <> __global__
void dev_coord_cart2pol<double>(const double* const arrayX,   const double* const arrayY,
                                      double* const arrayRad,       double* const arrayThe,
                                size_t length)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < length) {
        arrayRad[i] = sqrt(arrayX[i]*arrayX[i] + arrayY[i]*arrayY[i]);
        arrayThe[i] = atan2(arrayY[i],arrayX[i]);
    }
}

#ifndef NDEBUG
template <typename T> __global__
void dev_coord_cart2pol(const T* const arrayRow, const T* const arrayCol,
                              T* const arrayRad,       T* const arrayThe,
                        size_t length,
                        T      cRow, T      cCol,
                        size_t nRow, size_t nCol)
#else
template <typename T> __global__
void dev_coord_cart2pol(const T* const arrayRow, const T* const arrayCol,
                              T* const arrayRad,       T* const arrayThe,
                        size_t length,
                        T      cRow, T      cCol,
                        size_t , size_t )
#endif
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    T         x, y;

    if (i < length) {
        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(arrayRow[i] == ::round(arrayRow[i]));
        assert(arrayCol[i] == ::round(arrayCol[i]));

        x = arrayCol[i] - cCol;
        y = cRow - arrayRow[i];

        arrayRad[i] = sqrtf(x*x + y*y);
        arrayThe[i] = atan2f(y,x);
    }
}

#ifndef NDEBUG
template <> __global__
void dev_coord_cart2pol<double>(const double* const arrayRow, const double* const arrayCol,
                                      double* const arrayRad,       double* const arrayThe,
                                size_t length,
                                double cRow, double cCol,
                                size_t nRow, size_t nCol)
#else
template <> __global__
void dev_coord_cart2pol<double>(const double* const arrayRow, const double* const arrayCol,
                                      double* const arrayRad,       double* const arrayThe,
                                size_t length,
                                double cRow, double cCol,
                                size_t , size_t )
#endif

{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    double    x, y;

    if (i < length) {
        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(arrayRow[i] == ::round(arrayRow[i]));
        assert(arrayCol[i] == ::round(arrayCol[i]));

        x = arrayCol[i] - cCol;
        y = cRow - arrayRow[i];

        arrayRad[i] = sqrt(x*x + y*y);
        arrayThe[i] = atan2(y,x);
    }
}

template <typename T>
void cuda_coord_cart2pol(const T* const arrayX,   const T* const arrayY,
                               T* const arrayRad,       T* const arrayThe,
                         size_t length)
{
    assert(arrayX   != NULL && arrayY   != NULL);
    assert(arrayRad != NULL && arrayThe != NULL);
    assert(length > 0);

    dev_coord_cart2pol<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayX,   arrayY,
         arrayRad, arrayThe,
         length);
}

template <typename T>
void cuda_coord_cart2pol(const T* const arrayRow, const T* const arrayCol,
                               T* const arrayRad,       T* const arrayThe,
                         size_t length,
                         T      cRow, T      cCol,
                         size_t nRow, size_t nCol)
{
    assert(arrayRow != NULL && arrayCol != NULL);
    assert(arrayRad != NULL && arrayThe != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0);

    dev_coord_cart2pol<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayRow, arrayCol,
         arrayRad, arrayThe,
         length,
         cRow, cCol,
         nRow, nCol);
}

// instantiation
template
void cuda_coord_cart2pol<float >(const float*  const arrayRow, const float*  const arrayCol,
                                       float*  const arrayRad,       float*  const arrayThe,
                                 size_t length);
template
void cuda_coord_cart2pol<double>(const double* const arrayRow, const double* const arrayCol,
                                       double* const arrayRad,       double* const arrayThe,
                                 size_t length);
template
void cuda_coord_cart2pol<float >(const float*  const arrayRow, const float*  const arrayCol,
                                       float*  const arrayRad,       float*  const arrayThe,
                                 size_t length,
                                 float  cRow, float  cCol,
                                 size_t nRow, size_t nCol);
template
void cuda_coord_cart2pol<double>(const double* const arrayRow, const double* const arrayCol,
                                       double* const arrayRad,       double* const arrayThe,
                                 size_t length,
                                 double cRow, double cCol,
                                 size_t nRow, size_t nCol);

// pol2cart
template <typename T> __global__
void dev_coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                              T* const arrayX,         T* const arrayY,
                        size_t length)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < length) {
        assert(arrayRad[i] >= 0);

        arrayX[i] = arrayRad[i]*cosf(arrayThe[i]);
        arrayY[i] = arrayRad[i]*sinf(arrayThe[i]);
    }
}

template <> __global__
void dev_coord_pol2cart<double>(const double* const arrayRad, const double* const arrayThe,
                                      double* const arrayX,         double* const arrayY,
                                size_t length)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < length) {
        assert(arrayRad[i] >= 0);

        arrayX[i] = arrayRad[i]*cos(arrayThe[i]);
        arrayY[i] = arrayRad[i]*sin(arrayThe[i]);
    }
}

#ifndef NDEBUG
template <typename T> __global__
void dev_coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                              T* const arrayRow,       T* const arrayCol,
                        size_t length,
                        T      cRow, T      cCol,
                        size_t nRow, size_t nCol)
#else
template <typename T> __global__
void dev_coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                              T* const arrayRow,       T* const arrayCol,
                        size_t length,
                        T      cRow, T      cCol,
                        size_t , size_t )
#endif
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    T         tmpRow, tmpCol;

    if (i < length) {
        assert(arrayRad[i] >= 0);

        tmpRow = cRow - arrayRad[i]*sinf(arrayThe[i]);
        tmpCol = cCol + arrayRad[i]*cosf(arrayThe[i]);
        arrayRow[i] = rintf(tmpRow);
        arrayCol[i] = rintf(tmpCol);

        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(fabsf(arrayRow[i]-tmpRow) < 0.001);
        assert(fabsf(arrayCol[i]-tmpCol) < 0.001);
    }
}

#ifndef NDEBUG
template <> __global__
void dev_coord_pol2cart<double>(const double* const arrayRad, const double* const arrayThe,
                                      double* const arrayRow,       double* const arrayCol,
                                size_t length,
                                double cRow, double cCol,
                                size_t nRow, size_t nCol)
#else
template <> __global__
void dev_coord_pol2cart<double>(const double* const arrayRad, const double* const arrayThe,
                                      double* const arrayRow,       double* const arrayCol,
                                size_t length,
                                double cRow, double cCol,
                                size_t , size_t )
#endif
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    double    tmpRow, tmpCol;

    if (i < length) {
        assert(arrayRad[i] >= 0);

        tmpRow = cRow - arrayRad[i]*sin(arrayThe[i]);
        tmpCol = cCol + arrayRad[i]*cos(arrayThe[i]);
        arrayRow[i] = rint(tmpRow);
        arrayCol[i] = rint(tmpCol);

        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(fabs(arrayRow[i]-tmpRow) < 0.001);
        assert(fabs(arrayCol[i]-tmpCol) < 0.001);
    }
}

template <typename T>
void cuda_coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                               T* const arrayX,         T* const arrayY,
                         size_t length)
{
    assert(arrayRad != NULL && arrayThe != NULL);
    assert(arrayX   != NULL && arrayY   != NULL);
    assert(length > 0);

    dev_coord_pol2cart<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayRad, arrayThe,
         arrayX,   arrayY,
         length);
}

template <typename T>
void cuda_coord_pol2cart(const T* const arrayRad, const T* const arrayThe,
                               T* const arrayRow,       T* const arrayCol,
                         size_t length,
                         T      cRow, T      cCol,
                         size_t nRow, size_t nCol)
{
    assert(arrayRad != NULL && arrayThe != NULL);
    assert(arrayRow != NULL && arrayCol != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0);

    dev_coord_pol2cart<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayRad, arrayThe,
         arrayRow, arrayCol,
         length,
         cRow, cCol,
         nRow, nCol);
}

// instantiation
template
void cuda_coord_pol2cart<float >(const float*  const arrayRad, const float*  const arrayThe,
                                       float*  const arrayRow,       float*  const arrayCol,
                                 size_t length);
template
void cuda_coord_pol2cart<double>(const double* const arrayRad, const double* const arrayThe,
                                       double* const arrayRow,       double* const arrayCol,
                                 size_t length);
template
void cuda_coord_pol2cart<float >(const float*  const arrayRad, const float*  const arrayThe,
                                       float*  const arrayRow,       float*  const arrayCol,
                                 size_t length,
                                 float  cRow, float  cCol,
                                 size_t nRow, size_t nCol);
template
void cuda_coord_pol2cart<double>(const double* const arrayRad, const double* const arrayThe,
                                       double* const arrayRow,       double* const arrayCol,
                                 size_t length,
                                 double cRow, double cCol,
                                 size_t nRow, size_t nCol);

} // namespace gem

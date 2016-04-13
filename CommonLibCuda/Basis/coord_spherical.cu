/***********************************************************************
 *  File:       basis_spherical.cu
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
 * Cartesian - Spherical
 *************************/

// cart2sph
template <typename T> __global__
void dev_coord_cart2sph(const T* const arrayX,   const T* const arrayY,   const T* const arrayZ,
                              T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                        size_t length)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < length) {
        arrayRad[i] = sqrtf(arrayX[i]*arrayX[i] + arrayY[i]*arrayY[i] + arrayZ[i]*arrayZ[i]);
        arrayThe[i] = atan2f(sqrtf(arrayX[i]*arrayX[i] + arrayY[i]*arrayY[i]), arrayZ[i]);
        arrayPhi[i] = atan2f(arrayY[i],arrayX[i]);
    }
}

template <> __global__
void dev_coord_cart2sph<double>(const double* const arrayX,   const double* const arrayY,   const double* const arrayZ,
                                      double* const arrayRad,       double* const arrayThe,       double* const arrayPhi,
                                size_t length)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < length) {
        arrayRad[i] = sqrt(arrayX[i]*arrayX[i] + arrayY[i]*arrayY[i] + arrayZ[i]*arrayZ[i]);
        arrayThe[i] = atan2(sqrt(arrayX[i]*arrayX[i] + arrayY[i]*arrayY[i]), arrayZ[i]);
        arrayPhi[i] = atan2(arrayY[i],arrayX[i]);
    }
}

#ifndef NDEBUG
template <typename T> __global__
void dev_coord_cart2sph(const T* const arrayRow, const T* const arrayCol, const T* const arraySec,
                              T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                        size_t length,
                        T      cRow, T      cCol, T      cSec,
                        size_t nRow, size_t nCol, size_t nSec)
#else
template <typename T> __global__
void dev_coord_cart2sph(const T* const arrayRow, const T* const arrayCol, const T* const arraySec,
                              T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                        size_t length,
                        T      cRow, T      cCol, T      cSec,
                        size_t , size_t , size_t )
#endif
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    T         x, y, z;

    if (i < length) {
        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(arraySec[i] >= 0 && arraySec[i] < nSec);
        assert(arrayRow[i] == ::round(arrayRow[i]));
        assert(arrayCol[i] == ::round(arrayCol[i]));
        assert(arraySec[i] == ::round(arraySec[i]));

        x = cSec - arraySec[i];
        y = arrayCol[i] - cCol;
        z = cRow - arrayRow[i];

        arrayRad[i] = sqrtf(x*x + y*y + z*z);
        arrayThe[i] = atan2f(sqrtf(x*x + y*y), z);
        arrayPhi[i] = atan2f(y,x);
    }
}

#ifndef NDEBUG
template <> __global__
void dev_coord_cart2sph<double>(const double* const arrayRow, const double* const arrayCol, const double* const arraySec,
                                      double* const arrayRad,       double* const arrayThe,       double* const arrayPhi,
                                size_t length,
                                double cRow, double cCol, double cSec,
                                size_t nRow, size_t nCol, size_t nSec)
#else
template <> __global__
void dev_coord_cart2sph<double>(const double* const arrayRow, const double* const arrayCol, const double* const arraySec,
                                      double* const arrayRad,       double* const arrayThe,       double* const arrayPhi,
                                size_t length,
                                double cRow, double cCol, double cSec,
                                size_t , size_t , size_t )
#endif
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    double    x, y, z;

    if (i < length) {
        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(arraySec[i] >= 0 && arraySec[i] < nSec);
        assert(arrayRow[i] == ::round(arrayRow[i]));
        assert(arrayCol[i] == ::round(arrayCol[i]));
        assert(arraySec[i] == ::round(arraySec[i]));

        x = cSec - arraySec[i];
        y = arrayCol[i] - cCol;
        z = cRow - arrayRow[i];

        arrayRad[i] = sqrt(x*x + y*y + z*z);
        arrayThe[i] = atan2(sqrt(x*x + y*y), z);
        arrayPhi[i] = atan2(y,x);
    }
}

template <typename T>
void cuda_coord_cart2sph(const T* const arrayX,   const T* const arrayY,   const T* const arrayZ,
                               T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                         size_t length)
{
    assert(arrayX   != NULL && arrayY   != NULL && arrayZ   != NULL);
    assert(arrayRad != NULL && arrayThe != NULL && arrayPhi != NULL);
    assert(length > 0);

    dev_coord_cart2sph<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayX,   arrayY,   arrayZ,
         arrayRad, arrayThe, arrayPhi,
         length);
}

template <typename T>
void cuda_coord_cart2sph(const T* const arrayRow, const T* const arrayCol, const T* const arraySec,
                               T* const arrayRad,       T* const arrayThe,       T* const arrayPhi,
                         size_t length,
                         T      cRow, T      cCol, T      cSec,
                         size_t nRow, size_t nCol, size_t nSec)
{
    assert(arrayRow != NULL && arrayCol != NULL && arraySec != NULL);
    assert(arrayRad != NULL && arrayThe != NULL && arrayPhi != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    dev_coord_cart2sph<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayRow, arrayCol, arraySec,
         arrayRad, arrayThe, arrayPhi,
         length,
         cRow, cCol, cSec,
         nRow, nCol, nSec);
}

// instantiation
template
void cuda_coord_cart2sph<float >(const float*  const arrayRow, const float*  const arrayCol, const float*  const arraySec,
                                       float*  const arrayRad,       float*  const arrayThe,       float*  const arrayPhi,
                                 size_t length);
template
void cuda_coord_cart2sph<double>(const double* const arrayRow, const double* const arrayCol, const double* const arraySec,
                                       double* const arrayRad,       double* const arrayThe,       double* const arrayPhi,
                                 size_t length);
template
void cuda_coord_cart2sph<float >(const float*  const arrayRow, const float*  const arrayCol, const float*  const arraySec,
                                       float*  const arrayRad,       float*  const arrayThe,       float*  const arrayPhi,
                                 size_t length,
                                 float  cRow, float  cCol, float  cSec,
                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_coord_cart2sph<double>(const double* const arrayRow, const double* const arrayCol, const double* const arraySec,
                                       double* const arrayRad,       double* const arrayThe,       double* const arrayPhi,
                                 size_t length,
                                 double cRow, double cCol, double cSec,
                                 size_t nRow, size_t nCol, size_t nSec);

// sph2cart
template <typename T> __global__
void dev_coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                              T* const arrayX,         T* const arrayY,         T* const arrayZ,
                        size_t length)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    T         proj;

    if (i < length) {
        assert(arrayRad[i] >= 0);

        proj = arrayRad[i]*sinf(arrayThe[i]);

        arrayX[i] = proj*cosf(arrayPhi[i]);
        arrayY[i] = proj*sinf(arrayPhi[i]);
        arrayZ[i] = arrayRad[i]*cosf(arrayThe[i]);
    }
}

template <> __global__
void dev_coord_sph2cart<double>(const double* const arrayRad, const double* const arrayThe, const double* const arrayPhi,
                                      double* const arrayX,         double* const arrayY,         double* const arrayZ,
                                size_t length)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    double    proj;

    if (i < length) {
        assert(arrayRad[i] >= 0);

        proj = arrayRad[i]*sin(arrayThe[i]);

        arrayX[i] = proj*cos(arrayPhi[i]);
        arrayY[i] = proj*sin(arrayPhi[i]);
        arrayZ[i] = arrayRad[i]*cos(arrayThe[i]);
    }
}

#ifndef NDEBUG
template <typename T> __global__
void dev_coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                              T* const arrayRow,       T* const arrayCol,       T* const arraySec,
                        size_t length,
                        T      cRow, T      cCol, T      cSec,
                        size_t nRow, size_t nCol, size_t nSec)
#else
template <typename T> __global__
void dev_coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                              T* const arrayRow,       T* const arrayCol,       T* const arraySec,
                        size_t length,
                        T      cRow, T      cCol, T      cSec,
                        size_t , size_t , size_t )
#endif
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    T         proj, tmpRow, tmpCol, tmpSec;

    if (i < length) {
        assert(arrayRad[i] >= 0);

        proj = arrayRad[i]*sinf(arrayThe[i]);

        tmpRow = cRow - arrayRad[i]*cosf(arrayThe[i]);
        tmpCol = cCol + proj*sinf(arrayPhi[i]);
        tmpSec = cSec - proj*cosf(arrayPhi[i]);
        arrayRow[i] = rintf(tmpRow);
        arrayCol[i] = rintf(tmpCol);
        arraySec[i] = rintf(tmpSec);

        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(arraySec[i] >= 0 && arraySec[i] < nSec);
        assert(fabs(arrayRow[i]-tmpRow) < 0.001);
        assert(fabs(arrayCol[i]-tmpCol) < 0.001);
        assert(fabs(arraySec[i]-tmpSec) < 0.001);
    }
}

#ifndef NDEBUG
template <> __global__
void dev_coord_sph2cart<double>(const double* const arrayRad, const double* const arrayThe, const double* const arrayPhi,
                                      double* const arrayRow,       double* const arrayCol,       double* const arraySec,
                                size_t length,
                                double cRow, double cCol, double cSec,
                                size_t nRow, size_t nCol, size_t nSec)
#else
template <> __global__
void dev_coord_sph2cart<double>(const double* const arrayRad, const double* const arrayThe, const double* const arrayPhi,
                                      double* const arrayRow,       double* const arrayCol,       double* const arraySec,
                                size_t length,
                                double cRow, double cCol, double cSec,
                                size_t  , size_t , size_t )
#endif
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    double    proj, tmpRow, tmpCol, tmpSec;

    if (i < length) {
        assert(arrayRad[i] >= 0);

        proj = arrayRad[i]*sin(arrayThe[i]);

        tmpRow = cRow - arrayRad[i]*cos(arrayThe[i]);
        tmpCol = cCol + proj*sin(arrayPhi[i]);
        tmpSec = cSec - proj*cos(arrayPhi[i]);
        arrayRow[i] = rint(tmpRow);
        arrayCol[i] = rint(tmpCol);
        arraySec[i] = rint(tmpSec);

        assert(arrayRow[i] >= 0 && arrayRow[i] < nRow);
        assert(arrayCol[i] >= 0 && arrayCol[i] < nCol);
        assert(arraySec[i] >= 0 && arraySec[i] < nSec);
        assert(fabs(arrayRow[i]-tmpRow) < 0.001);
        assert(fabs(arrayCol[i]-tmpCol) < 0.001);
        assert(fabs(arraySec[i]-tmpSec) < 0.001);
    }
}

template <typename T>
void cuda_coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                               T* const arrayX,         T* const arrayY,         T* const arrayZ,
                         size_t length)
{
    assert(arrayRad != NULL && arrayThe != NULL && arrayPhi != NULL);
    assert(arrayX   != NULL && arrayY   != NULL && arrayZ   != NULL);
    assert(length > 0);

    dev_coord_sph2cart<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayRad, arrayThe, arrayPhi,
         arrayX,   arrayY,   arrayZ,
         length);
}

template <typename T>
void cuda_coord_sph2cart(const T* const arrayRad, const T* const arrayThe, const T* const arrayPhi,
                               T* const arrayRow,       T* const arrayCol,       T* const arraySec,
                         size_t length,
                         T      cRow, T      cCol, T      cSec,
                         size_t nRow, size_t nCol, size_t nSec)
{
    assert(arrayRad != NULL && arrayThe != NULL && arrayPhi != NULL);
    assert(arrayRow != NULL && arrayCol != NULL && arraySec != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    dev_coord_sph2cart<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayRad, arrayThe, arrayPhi,
         arrayRow, arrayCol, arraySec,
         length,
         cRow, cCol, cSec,
         nRow, nCol, nSec);
}

// instantiation
template
void cuda_coord_sph2cart<float >(const float*  const arrayRad, const float*  const arrayThe, const float*  const arrayPhi,
                                       float*  const arrayRow,       float*  const arrayCol,       float*  const arraySec,
                                 size_t length);
template
void cuda_coord_sph2cart<double>(const double* const arrayRad, const double* const arrayThe, const double* const arrayPhi,
                                       double* const arrayRow,       double* const arrayCol,       double* const arraySec,
                                 size_t length);
template
void cuda_coord_sph2cart<float >(const float*  const arrayRad, const float*  const arrayThe, const float*  const arrayPhi,
                                       float*  const arrayRow,       float*  const arrayCol,       float*  const arraySec,
                                 size_t length,
                                 float  cRow, float  cCol, float  cSec,
                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_coord_sph2cart<double>(const double* const arrayRad, const double* const arrayThe, const double* const arrayPhi,
                                       double* const arrayRow,       double* const arrayCol,       double* const arraySec,
                                 size_t length,
                                 double cRow, double cCol, double cSec,
                                 size_t nRow, size_t nCol, size_t nSec);

} // namespace gem

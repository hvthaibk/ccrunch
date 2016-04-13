/***********************************************************************
 *  File:       transform_rottrans.cu
 *
 *  Purpose:    Implementation of transformation functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "transform.cuh"

namespace gem {

/**************************
 * Rotation-Translation
 *************************/

// 2D
template <typename T> __global__
void dev_transform_rotate_translate(const T* const arraySrc, T* const arrayDst,
                                    size_t nRow, size_t nCol,
                                    T      cRow, T      cCol,
                                    T    cTheta, T    sTheta,
                                    T nRowOff,   T nColOff)
{
    size_t    iRow = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iCol = blockDim.x * blockIdx.x + threadIdx.x;

    ptrdiff_t    intRow, intCol, intRow1, intCol1;
    bool         bRow, bCol, bRow1, bCol1;
    T            iRowOri, iColOri, iRowRot, iColRot;
    T            v00, v01, v10, v11;
    T            fracRow, fracCol;

    if (iRow < nRow && iCol < nCol) {
        iRowOri = (T) iRow - cRow - nRowOff;
        iColOri = (T) iCol - cCol - nColOff;

        iRowRot =  cTheta*iRowOri + sTheta*iColOri + cRow;
        iColRot = -sTheta*iRowOri + cTheta*iColOri + cCol;

        intRow = (ptrdiff_t) std::floor(iRowRot);
        intCol = (ptrdiff_t) std::floor(iColRot);

        intRow1 = intRow + 1;
        intCol1 = intCol + 1;

        bRow  = (intRow  >= 0) && (intRow  <= (ptrdiff_t) nRow-1);
        bCol  = (intCol  >= 0) && (intCol  <= (ptrdiff_t) nCol-1);
        bRow1 = (intRow1 >= 0) && (intRow1 <= (ptrdiff_t) nRow-1);
        bCol1 = (intCol1 >= 0) && (intCol1 <= (ptrdiff_t) nCol-1);

        v00 = 0;    v01 = 0;    v10 = 0;    v11 = 0;

        if (bRow  && bCol ) v00 = arraySrc[intRow *nCol+intCol ];
        if (bRow  && bCol1) v01 = arraySrc[intRow *nCol+intCol1];
        if (bRow1 && bCol ) v10 = arraySrc[intRow1*nCol+intCol ];
        if (bRow1 && bCol1) v11 = arraySrc[intRow1*nCol+intCol1];

        fracRow = iRowRot - (T) intRow;
        fracCol = iColRot - (T) intCol;

        arrayDst[iRow*nCol+iCol] = (1-fracRow) * (1-fracCol) * v00 +
                                   (1-fracRow) * fracCol     * v01 +
                                   fracRow     * (1-fracCol) * v10 +
                                   fracRow     * fracCol     * v11;
    }
}

template <> __global__
void dev_transform_rotate_translate<double>(const double* const arraySrc, double* const arrayDst,
                                            size_t nRow,    size_t nCol,
                                            double cRow,    double cCol,
                                            double cTheta,  double sTheta,
                                            double nRowOff, double nColOff)
{
    size_t    iRow = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iCol = blockDim.x * blockIdx.x + threadIdx.x;

    ptrdiff_t    intRow, intCol, intRow1, intCol1;
    bool         bRow, bCol, bRow1, bCol1;
    double       iRowOri, iColOri, iRowRot, iColRot;
    double       v00, v01, v10, v11;
    double       fracRow, fracCol;

    if (iRow < nRow && iCol < nCol) {
        iRowOri = (double) iRow - cRow - nRowOff;
        iColOri = (double) iCol - cCol - nColOff;

        iRowRot =  cTheta*iRowOri + sTheta*iColOri + cRow;
        iColRot = -sTheta*iRowOri + cTheta*iColOri + cCol;

        intRow = (ptrdiff_t) floor(iRowRot);
        intCol = (ptrdiff_t) floor(iColRot);

        intRow1 = intRow + 1;
        intCol1 = intCol + 1;

        bRow  = (intRow  >= 0) && (intRow  <= (ptrdiff_t) nRow-1);
        bCol  = (intCol  >= 0) && (intCol  <= (ptrdiff_t) nCol-1);
        bRow1 = (intRow1 >= 0) && (intRow1 <= (ptrdiff_t) nRow-1);
        bCol1 = (intCol1 >= 0) && (intCol1 <= (ptrdiff_t) nCol-1);

        v00 = 0;    v01 = 0;    v10 = 0;    v11 = 0;

        if (bRow  && bCol ) v00 = arraySrc[intRow *nCol+intCol ];
        if (bRow  && bCol1) v01 = arraySrc[intRow *nCol+intCol1];
        if (bRow1 && bCol ) v10 = arraySrc[intRow1*nCol+intCol ];
        if (bRow1 && bCol1) v11 = arraySrc[intRow1*nCol+intCol1];

        fracRow = iRowRot - (double) intRow;
        fracCol = iColRot - (double) intCol;

        arrayDst[iRow*nCol+iCol] = (1-fracRow) * (1-fracCol) * v00 +
                                   (1-fracRow) * fracCol     * v01 +
                                   fracRow     * (1-fracCol) * v10 +
                                   fracRow     * fracCol     * v11;
    }
}

template <typename T>
void cuda_transform_rotate_translate(const T* const arraySrc, T* const arrayDst,
                                     size_t nRow,        size_t nCol,
                                     T      theta,
                                     T      nRowOff,     T      nColOff,
                                     eInter inter)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);

    switch (inter) {
        case INTER_NEAREST:
            ERROR("cuda_transform_rotate", "unsupported interpolation mode");
            break;
        case INTER_LINEAR:
            break;
        case INTER_CUBIC:
            ERROR("cuda_transform_rotate", "unsupported interpolation mode");
            break;
        default:
            ERROR("cuda_transform_rotate", "unsupported interpolation mode");
    }

    T       cTheta = std::cos(theta);
    T       sTheta = std::sin(theta);
    T       cRow = transform_centerCoord<T>(nRow);
    T       cCol = transform_centerCoord<T>(nCol);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nCol, dimBlock.x),
                     iDivUp(nRow, dimBlock.y));

    dev_transform_rotate_translate<<<dimGrid, dimBlock>>>
        (arraySrc, arrayDst,
         nRow, nCol,
         cRow, cCol,
         cTheta,  sTheta,
         nRowOff, nColOff);
}

// instantiation
template
void cuda_transform_rotate_translate<float >(const float*  const arraySrc, float*  const arrayDst,
                                             size_t nRow,        size_t nCol,
                                             float  theta,
                                             float  nRowOff,     float  nColOff,
                                             eInter inter);
template
void cuda_transform_rotate_translate<double>(const double* const arraySrc, double* const arrayDst,
                                             size_t nRow,        size_t nCol,
                                             double theta,
                                             double nRowOff,     double nColOff,
                                             eInter inter);

// 3D
template <typename T> __global__
void dev_transform_rotate_translate(const T* const arraySrc, T* const arrayDst,
                                    size_t nRow, size_t nCol, size_t nSec,
                                    T      cRow, T      cCol, T      cSec,
                                    T m11, T m12, T m13,
                                    T m21, T m22, T m23,
                                    T m31, T m32, T m33,
                                    T nRowOff, T nColOff, T nSecOff)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#else
    size_t    iRow = blockDim.z * blockIdx.z + threadIdx.z;
    size_t    iCol = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSec = blockDim.x * blockIdx.x + threadIdx.x;

    ptrdiff_t    intRow, intCol, intSec, intRow1, intCol1, intSec1;
    bool         bRow, bCol, bSec, bRow1, bCol1, bSec1;
    T            iRowOri, iColOri, iSecOri, iRowRot, iColRot, iSecRot;
    T            v000, v001, v010, v011, v100, v101, v110, v111;
    T            fracRow, fracCol, fracSec;

    if (iRow < nRow && iCol < nCol && iSec <nSec) {
        iRowOri = (T) iRow - cRow - nRowOff;
        iColOri = (T) iCol - cCol - nColOff;
        iSecOri = (T) iSec - cSec - nSecOff;

        iRowRot = m11*iRowOri + m12*iColOri + m13*iSecOri + cRow;
        iColRot = m21*iRowOri + m22*iColOri + m23*iSecOri + cCol;
        iSecRot = m31*iRowOri + m32*iColOri + m33*iSecOri + cSec;

        intRow = (ptrdiff_t) std::floor(iRowRot);
        intCol = (ptrdiff_t) std::floor(iColRot);
        intSec = (ptrdiff_t) std::floor(iSecRot);

        intRow1 = intRow + 1;
        intCol1 = intCol + 1;
        intSec1 = intSec + 1;

        bRow  = (intRow  >= 0) && (intRow  <= (ptrdiff_t) nRow-1);
        bCol  = (intCol  >= 0) && (intCol  <= (ptrdiff_t) nCol-1);
        bSec  = (intSec  >= 0) && (intSec  <= (ptrdiff_t) nSec-1);
        bRow1 = (intRow1 >= 0) && (intRow1 <= (ptrdiff_t) nRow-1);
        bCol1 = (intCol1 >= 0) && (intCol1 <= (ptrdiff_t) nCol-1);
        bSec1 = (intSec1 >= 0) && (intSec1 <= (ptrdiff_t) nSec-1);

        v000 = 0;    v001 = 0;    v010 = 0;    v011 = 0;
        v100 = 0;    v101 = 0;    v110 = 0;    v111 = 0;

        if (bRow  && bCol  && bSec ) v000 = arraySrc[intRow *nCol*nSec+intCol *nSec+intSec ];
        if (bRow  && bCol  && bSec1) v001 = arraySrc[intRow *nCol*nSec+intCol *nSec+intSec1];
        if (bRow  && bCol1 && bSec ) v010 = arraySrc[intRow *nCol*nSec+intCol1*nSec+intSec ];
        if (bRow  && bCol1 && bSec1) v011 = arraySrc[intRow *nCol*nSec+intCol1*nSec+intSec1];
        if (bRow1 && bCol  && bSec ) v100 = arraySrc[intRow1*nCol*nSec+intCol *nSec+intSec ];
        if (bRow1 && bCol  && bSec1) v101 = arraySrc[intRow1*nCol*nSec+intCol *nSec+intSec1];
        if (bRow1 && bCol1 && bSec ) v110 = arraySrc[intRow1*nCol*nSec+intCol1*nSec+intSec ];
        if (bRow1 && bCol1 && bSec1) v111 = arraySrc[intRow1*nCol*nSec+intCol1*nSec+intSec1];

        fracRow = iRowRot - (T) intRow;
        fracCol = iColRot - (T) intCol;
        fracSec = iSecRot - (T) intSec;

        arrayDst[(iRow*nCol+iCol)*nSec+iSec]
                    = (1-fracRow) * (1-fracCol) * (1-fracSec) * v000
                    + (1-fracRow) * (1-fracCol) * fracSec     * v001
                    + (1-fracRow) * fracCol     * (1-fracSec) * v010
                    + (1-fracRow) * fracCol     * fracSec     * v011
                    + fracRow     * (1-fracCol) * (1-fracSec) * v100
                    + fracRow     * (1-fracCol) * fracSec     * v101
                    + fracRow     * fracCol     * (1-fracSec) * v110
                    + fracRow     * fracCol     * fracSec     * v111;
    }
#endif
}

template <> __global__
void dev_transform_rotate_translate<double>(const double* const arraySrc, double* const arrayDst,
                                            size_t nRow, size_t nCol, size_t nSec,
                                            double cRow, double cCol, double cSec,
                                            double m11, double m12, double m13,
                                            double m21, double m22, double m23,
                                            double m31, double m32, double m33,
                                            double nRowOff, double nColOff, double nSecOff)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#else
    size_t    iRow = blockDim.z * blockIdx.z + threadIdx.z;
    size_t    iCol = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSec = blockDim.x * blockIdx.x + threadIdx.x;

    ptrdiff_t    intRow, intCol, intSec, intRow1, intCol1, intSec1;
    bool         bRow, bCol, bSec, bRow1, bCol1, bSec1;
    double       iRowOri, iColOri, iSecOri, iRowRot, iColRot, iSecRot;
    double       v000, v001, v010, v011, v100, v101, v110, v111;
    double       fracRow, fracCol, fracSec;

    if (iRow < nRow && iCol < nCol && iSec <nSec) {
        iRowOri = (double) iRow - cRow - nRowOff;
        iColOri = (double) iCol - cCol - nColOff;
        iSecOri = (double) iSec - cSec - nSecOff;

        iRowRot = m11*iRowOri + m12*iColOri + m13*iSecOri + cRow;
        iColRot = m21*iRowOri + m22*iColOri + m23*iSecOri + cCol;
        iSecRot = m31*iRowOri + m32*iColOri + m33*iSecOri + cSec;

        intRow = (ptrdiff_t) floor(iRowRot);
        intCol = (ptrdiff_t) floor(iColRot);
        intSec = (ptrdiff_t) floor(iSecRot);

        intRow1 = intRow + 1;
        intCol1 = intCol + 1;
        intSec1 = intSec + 1;

        bRow  = (intRow  >= 0) && (intRow  <= (ptrdiff_t) nRow-1);
        bCol  = (intCol  >= 0) && (intCol  <= (ptrdiff_t) nCol-1);
        bSec  = (intSec  >= 0) && (intSec  <= (ptrdiff_t) nSec-1);
        bRow1 = (intRow1 >= 0) && (intRow1 <= (ptrdiff_t) nRow-1);
        bCol1 = (intCol1 >= 0) && (intCol1 <= (ptrdiff_t) nCol-1);
        bSec1 = (intSec1 >= 0) && (intSec1 <= (ptrdiff_t) nSec-1);

        v000 = 0;    v001 = 0;    v010 = 0;    v011 = 0;
        v100 = 0;    v101 = 0;    v110 = 0;    v111 = 0;

        if (bRow  && bCol  && bSec ) v000 = arraySrc[intRow *nCol*nSec+intCol *nSec+intSec ];
        if (bRow  && bCol  && bSec1) v001 = arraySrc[intRow *nCol*nSec+intCol *nSec+intSec1];
        if (bRow  && bCol1 && bSec ) v010 = arraySrc[intRow *nCol*nSec+intCol1*nSec+intSec ];
        if (bRow  && bCol1 && bSec1) v011 = arraySrc[intRow *nCol*nSec+intCol1*nSec+intSec1];
        if (bRow1 && bCol  && bSec ) v100 = arraySrc[intRow1*nCol*nSec+intCol *nSec+intSec ];
        if (bRow1 && bCol  && bSec1) v101 = arraySrc[intRow1*nCol*nSec+intCol *nSec+intSec1];
        if (bRow1 && bCol1 && bSec ) v110 = arraySrc[intRow1*nCol*nSec+intCol1*nSec+intSec ];
        if (bRow1 && bCol1 && bSec1) v111 = arraySrc[intRow1*nCol*nSec+intCol1*nSec+intSec1];

        fracRow = iRowRot - (double) intRow;
        fracCol = iColRot - (double) intCol;
        fracSec = iSecRot - (double) intSec;

        arrayDst[(iRow*nCol+iCol)*nSec+iSec]
                    = (1-fracRow) * (1-fracCol) * (1-fracSec) * v000
                    + (1-fracRow) * (1-fracCol) * fracSec     * v001
                    + (1-fracRow) * fracCol     * (1-fracSec) * v010
                    + (1-fracRow) * fracCol     * fracSec     * v011
                    + fracRow     * (1-fracCol) * (1-fracSec) * v100
                    + fracRow     * (1-fracCol) * fracSec     * v101
                    + fracRow     * fracCol     * (1-fracSec) * v110
                    + fracRow     * fracCol     * fracSec     * v111;
    }
#endif
}

template <typename T>
void cuda_transform_rotate_translate(const T* const arraySrc, T* const arrayDst,
                                     size_t nRow,    size_t nCol,    size_t nSec,
                                     T      alpha,   T      beta,    T      gamma,
                                     T      nRowOff, T      nColOff, T      nSecOff,
                                     eInter inter)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    switch (inter) {
        case INTER_NEAREST:
            ERROR("cuda_transform_rotate_translate", "unsupported interpolation mode");
            break;
        case INTER_LINEAR:
            break;
        case INTER_CUBIC:
            ERROR("cuda_transform_rotate_translate", "unsupported interpolation mode");
            break;
        default:
            ERROR("cuda_transform_rotate_translate", "unsupported interpolation mode");
    }

    T       m11, m12, m13, m21, m22, m23, m31, m32, m33;
    T       cRow, cCol, cSec;

    cRow = transform_centerCoord<T>(nRow);
    cCol = transform_centerCoord<T>(nCol);
    cSec = transform_centerCoord<T>(nSec);

    transform_rotmat(-alpha, -beta, -gamma,
                     m11, m12, m13, m21, m22, m23, m31, m32, m33,
                     ROT3D_RIGHT_ZYZ);

#ifdef __GEM_CUDA_ARCH_HOST_130__
    ERROR("cuda_transform_rotate_translate", "unsupported hardware");
#else
    dim3    dimBlock(BLOCK_3D_NSEC, BLOCK_3D_NCOL, BLOCK_3D_NROW);
    dim3    dimGrid (iDivUp(nSec, dimBlock.x),
                     iDivUp(nCol, dimBlock.y),
                     iDivUp(nRow, dimBlock.z));
#endif

    dev_transform_rotate_translate<<<dimGrid, dimBlock>>>
        (arraySrc, arrayDst,
         nRow, nCol, nSec,
         cRow, cCol, cSec,
         m11, m12, m13, m21, m22, m23, m31, m32, m33,
         nRowOff, nColOff, nSecOff);
}

// instantiation
template
void cuda_transform_rotate_translate<float >(const float*  const arraySrc, float * const arrayDst,
                                             size_t nRow,    size_t nCol,    size_t nSec,
                                             float  alpha,   float  beta,    float  gamma,
                                             float  nRowOff, float  nColOff, float  nSecOff,
                                             eInter inter);
template
void cuda_transform_rotate_translate<double>(const double* const arraySrc, double* const arrayDst,
                                             size_t nRow,    size_t nCol,    size_t nSec,
                                             double alpha,   double beta,    double gamma,
                                             double nRowOff, double nColOff, double nSecOff,
                                             eInter inter);

} // namespace gem

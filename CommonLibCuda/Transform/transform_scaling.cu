/***********************************************************************
 *  File:       transform_scaling.cu
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
 * Scaling
 *************************/

// 1D
template <typename T> __global__
void dev_transform_scale(const T* const arraySrc, T* const arrayDst,
                         size_t nRow,
                         T      cRow,
                         size_t nRowDst,
                         T      cRowDst,
                         T factor)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    ptrdiff_t    intRow, intRow1;
    bool         bRow, bRow1;
    T            iRowScale;
    T            v0, v1;
    T            fracRow;

    if (iRow < nRowDst) {
        iRowScale = (iRow - cRowDst) / factor + cRow;

        intRow = (ptrdiff_t) std::floor(iRowScale);

        intRow1 = intRow + 1;

        bRow  = (intRow  >= 0) && (intRow  <= (ptrdiff_t) nRow-1);
        bRow1 = (intRow1 >= 0) && (intRow1 <= (ptrdiff_t) nRow-1);

        v0 = 0;     v1 = 0;

        if (bRow )  v0 = arraySrc[intRow ];
        if (bRow1)  v1 = arraySrc[intRow1];

        fracRow = iRowScale - (T) intRow;

        arrayDst[iRow] = (1-fracRow) * v0 + fracRow * v1;
    }
}

template <typename T>
void cuda_transform_scale(const T* const arraySrc, T* const arrayDst,
                          size_t nRow,
                          T factor,
                          eInter inter)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0 );

    switch (inter) {
        case INTER_NEAREST:
            ERROR("cuda_transform_scale", "unsupported interpolation mode");
            break;
        case INTER_LINEAR:
            break;
        case INTER_CUBIC:
            ERROR("cuda_transform_scale", "unsupported interpolation mode");
            break;
        default:
            ERROR("cuda_transform_scale", "unsupported interpolation mode");
    }

    T            cRow;
    size_t       nRowDst;
    T            cRowDst;

    cRow = transform_centerCoord<T>(nRow);

    nRowDst = transform_scaleSize(nRow, factor);

    cRowDst = transform_centerCoord<T>(nRowDst);

    dev_transform_scale<<<iDivUp(nRowDst, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arraySrc, arrayDst,
         nRow,
         cRow,
         nRowDst,
         cRowDst,
         factor);
}

// instantiation
template
void cuda_transform_scale<float >(const float*  const arraySrc, float * const arrayDst,
                                  size_t nRow,
                                  float  factor,
                                  eInter inter);
template
void cuda_transform_scale<double>(const double* const arraySrc, double* const arrayDst,
                                  size_t nRow,
                                  double factor,
                                  eInter inter);

// 2D
template <typename T> __global__
void dev_transform_scale(const T* const arraySrc, T* const arrayDst,
                         size_t nRow,    size_t nCol,
                         T      cRow,    T      cCol,
                         size_t nRowDst, size_t nColDst,
                         T      cRowDst, T      cColDst,
                         T factor)
{
    size_t    iRow = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iCol = blockDim.x * blockIdx.x + threadIdx.x;

    ptrdiff_t    intRow, intCol, intRow1, intCol1;
    bool         bRow, bCol, bRow1, bCol1;
    T            iRowScale, iColScale;
    T            v00, v01, v10, v11;
    T            fracRow, fracCol;

    if (iRow < nRowDst && iCol < nColDst) {
        iRowScale = (iRow - cRowDst) / factor + cRow;
        iColScale = (iCol - cColDst) / factor + cCol;

        intRow = (ptrdiff_t) std::floor(iRowScale);
        intCol = (ptrdiff_t) std::floor(iColScale);

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

        fracRow = iRowScale - (T) intRow;
        fracCol = iColScale - (T) intCol;

        arrayDst[iRow*nColDst+iCol] = (1-fracRow) * (1-fracCol) * v00 +
                                      (1-fracRow) * fracCol     * v01 +
                                      fracRow     * (1-fracCol) * v10 +
                                      fracRow     * fracCol     * v11;
    }
}

template <typename T>
void cuda_transform_scale(const T* const arraySrc, T* const arrayDst,
                          size_t nRow, size_t nCol,
                          T factor,
                          eInter inter)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);

    switch (inter) {
        case INTER_NEAREST:
            ERROR("cuda_transform_scale", "unsupported interpolation mode");
            break;
        case INTER_LINEAR:
            break;
        case INTER_CUBIC:
            ERROR("cuda_transform_scale", "unsupported interpolation mode");
            break;
        default:
            ERROR("cuda_transform_scale", "unsupported interpolation mode");
    }

    T            cRow, cCol;
    size_t       nRowDst, nColDst;
    T            cRowDst, cColDst;

    cRow = transform_centerCoord<T>(nRow);
    cCol = transform_centerCoord<T>(nCol);

    nRowDst = transform_scaleSize(nRow, factor);
    nColDst = transform_scaleSize(nCol, factor);

    cRowDst = transform_centerCoord<T>(nRowDst);
    cColDst = transform_centerCoord<T>(nColDst);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nColDst, dimBlock.x),
                     iDivUp(nRowDst, dimBlock.y));

    dev_transform_scale<<<dimGrid, dimBlock>>>
        (arraySrc, arrayDst,
         nRow,    nCol,
         cRow,    cCol,
         nRowDst, nColDst,
         cRowDst, cColDst,
         factor);
}

// instantiation
template
void cuda_transform_scale<float >(const float*  const arraySrc, float * const arrayDst,
                                  size_t nRow, size_t nCol,
                                  float  factor,
                                  eInter inter);
template
void cuda_transform_scale<double>(const double* const arraySrc, double* const arrayDst,
                                  size_t nRow, size_t nCol,
                                  double factor,
                                  eInter inter);

// 3D
template <typename T> __global__
void dev_transform_scale(const T* const arraySrc, T* const arrayDst,
                         size_t nRow,    size_t nCol,    size_t nSec,
                         T      cRow,    T      cCol,    T      cSec,
                         size_t nRowDst, size_t nColDst, size_t nSecDst,
                         T      cRowDst, T      cColDst, T      cSecDst,
                         T factor)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#else
    size_t    iRow = blockDim.z * blockIdx.z + threadIdx.z;
    size_t    iCol = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSec = blockDim.x * blockIdx.x + threadIdx.x;

    ptrdiff_t    intRow, intCol, intSec, intRow1, intCol1, intSec1;
    bool         bRow, bCol, bSec, bRow1, bCol1, bSec1;
    T            iRowScale, iColScale, iSecScale;
    T            v000, v001, v010, v011, v100, v101, v110, v111;
    T            fracRow, fracCol, fracSec;

    if (iRow < nRowDst && iCol < nColDst && iSec <nSecDst) {
        iRowScale = (iRow - cRowDst) / factor + cRow;
        iColScale = (iCol - cColDst) / factor + cCol;
        iSecScale = (iSec - cSecDst) / factor + cSec;

        intRow = (ptrdiff_t) std::floor(iRowScale);
        intCol = (ptrdiff_t) std::floor(iColScale);
        intSec = (ptrdiff_t) std::floor(iSecScale);

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

        fracRow = iRowScale - (T) intRow;
        fracCol = iColScale - (T) intCol;
        fracSec = iSecScale - (T) intSec;

        arrayDst[(iRow*nColDst+iCol)*nSecDst+iSec]
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
void cuda_transform_scale(const T* const arraySrc, T* const arrayDst,
                          size_t nRow, size_t nCol, size_t nSec,
                          T factor,
                          eInter inter)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    switch (inter) {
        case INTER_NEAREST:
            ERROR("cuda_transform_scale", "unsupported interpolation mode");
            break;
        case INTER_LINEAR:
            break;
        case INTER_CUBIC:
            ERROR("cuda_transform_scale", "unsupported interpolation mode");
            break;
        default:
            ERROR("cuda_transform_scale", "unsupported interpolation mode");
    }

    T            cRow, cCol, cSec;
    size_t       nRowDst, nColDst, nSecDst;
    T            cRowDst, cColDst, cSecDst;

    cRow = transform_centerCoord<T>(nRow);
    cCol = transform_centerCoord<T>(nCol);
    cSec = transform_centerCoord<T>(nSec);

    nRowDst = transform_scaleSize(nRow, factor);
    nColDst = transform_scaleSize(nCol, factor);
    nSecDst = transform_scaleSize(nSec, factor);

    cRowDst = transform_centerCoord<T>(nRowDst);
    cColDst = transform_centerCoord<T>(nColDst);
    cSecDst = transform_centerCoord<T>(nSecDst);

#ifdef __GEM_CUDA_ARCH_HOST_130__
    ERROR("cuda_transform_scale", "unsupported hardware");
    /*dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nSecDst, dimBlock.x),
                     iDivUp(nColDst, dimBlock.y));*/
#else
    dim3    dimBlock(BLOCK_3D_NSEC, BLOCK_3D_NCOL, BLOCK_3D_NROW);
    dim3    dimGrid (iDivUp(nSecDst, dimBlock.x),
                     iDivUp(nColDst, dimBlock.y),
                     iDivUp(nRowDst, dimBlock.z));
#endif

    dev_transform_scale<<<dimGrid, dimBlock>>>
        (arraySrc, arrayDst,
         nRow,    nCol,    nSec,
         cRow,    cCol,    cSec,
         nRowDst, nColDst, nSecDst,
         cRowDst, cColDst, cSecDst,
         factor);
}

// instantiation
template
void cuda_transform_scale<float >(const float*  const arraySrc, float * const arrayDst,
                                  size_t nRow, size_t nCol, size_t nSec,
                                  float  factor,
                                  eInter inter);
template
void cuda_transform_scale<double>(const double* const arraySrc, double* const arrayDst,
                                  size_t nRow, size_t nCol, size_t nSec,
                                  double factor,
                                  eInter inter);

} // namespace gem

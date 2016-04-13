/***********************************************************************
 *  File:       array_shifting.cu
 *
 *  Purpose:    Implementation of array-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "array.cuh"

namespace gem {

/*****************************************
 * Shifting
 ****************************************/

// 1D
template <typename T> __global__
void dev_array_shift(const T* const arraySrc, T* const arrayDst,
                     size_t    nRow,
                     ptrdiff_t nRowOff)
{
    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    iRowDst = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iRowSrc;

    if (iRowDst < nRowSigned) {
        iRowSrc = iRowDst - nRowOff;

        if ((iRowSrc >= 0) && (iRowSrc < nRowSigned)) {
            arrayDst[iRowDst] = arraySrc[iRowSrc];
        }
        else {
            arrayDst[iRowDst] = 0;
        }
    }
}

template <typename T>
void cuda_array_shift(const T* const arraySrc, T* const arrayDst,
                      size_t    nRow,
                      ptrdiff_t nRowOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_shift<<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arraySrc, arrayDst, nRow, nRowOff);
}

// instantiation
template
void cuda_array_shift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                                size_t    nRow,
                                ptrdiff_t nRowOff);
template
void cuda_array_shift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                                size_t    nRow,
                                ptrdiff_t nRowOff);
template
void cuda_array_shift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                                size_t    nRow,
                                ptrdiff_t nRowOff);
template
void cuda_array_shift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                                size_t    nRow,
                                ptrdiff_t nRowOff);
template
void cuda_array_shift<float   >(const float*    const arraySrc, float*    const arrayDst,
                                size_t    nRow,
                                ptrdiff_t nRowOff);
template
void cuda_array_shift<double  >(const double*   const arraySrc, double*   const arrayDst,
                                size_t    nRow,
                                ptrdiff_t nRowOff);

// 2D
template <typename T> __global__
void dev_array_shift(const T* const arraySrc, T* const arrayDst,
                     size_t    nRow,    size_t    nCol,
                     ptrdiff_t nRowOff, ptrdiff_t nColOff)
{
    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    nColSigned = (ptrdiff_t) nCol;
    ptrdiff_t    iRowDst = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iColDst = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iRowSrc, iColSrc;

    if (iRowDst < nRowSigned && iColDst < nColSigned) {
        iRowSrc = iRowDst - nRowOff;
        iColSrc = iColDst - nColOff;

        if ((iRowSrc >= 0) && (iRowSrc < nRowSigned) &&
            (iColSrc >= 0) && (iColSrc < nColSigned)) {
            arrayDst[iRowDst*nColSigned+iColDst] =
            arraySrc[iRowSrc*nColSigned+iColSrc];
        }
        else {
            arrayDst[iRowDst*nColSigned+iColDst] = 0;
        }
    }
}

template <typename T>
void cuda_array_shift(const T* const arraySrc, T* const arrayDst,
                      size_t    nRow,    size_t    nCol,
                      ptrdiff_t nRowOff, ptrdiff_t nColOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nCol, dimBlock.x),
                     iDivUp(nRow, dimBlock.y));

    dev_array_shift<<<dimGrid, dimBlock>>>
        (arraySrc, arrayDst, nRow, nCol, nRowOff, nColOff);
}

// instantiation
template
void cuda_array_shift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                                size_t    nRow,    size_t    nCol,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_shift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                                size_t    nRow,    size_t    nCol,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_shift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                                size_t    nRow,    size_t    nCol,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_shift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                                size_t    nRow,    size_t    nCol,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_shift<float   >(const float*    const arraySrc, float*    const arrayDst,
                                size_t    nRow,    size_t    nCol,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_shift<double  >(const double*   const arraySrc, double*   const arrayDst,
                                size_t    nRow,    size_t    nCol,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D
template <typename T> __global__
void dev_array_shift(const T* const arraySrc, T* const arrayDst,
                     size_t    nRow,    size_t    nCol,    size_t    nSec,
                     ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    nColSigned = (ptrdiff_t) nCol;
    ptrdiff_t    nSecSigned = (ptrdiff_t) nSec;
    ptrdiff_t    iRowDst;
    ptrdiff_t    iColDst = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iSecDst = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iRowSrc, iColSrc, iSecSrc;

    if (iColDst < nColSigned & iSecDst < nSecSigned) {
        iColSrc = iColDst - nColOff;
        iSecSrc = iSecDst - nSecOff;

        for (iRowDst = 0; iRowDst < nRowSigned; iRowDst++) {
            iRowSrc = iRowDst - nRowOff;

            if ((iRowSrc >= 0) && (iRowSrc < nRowSigned) &&
                (iColSrc >= 0) && (iColSrc < nColSigned) &&
                (iSecSrc >= 0) && (iSecSrc < nSecSigned)) {
                arrayDst[(iRowDst*nColSigned+iColDst)*nSecSigned+iSecDst] =
                arraySrc[(iRowSrc*nColSigned+iColSrc)*nSecSigned+iSecSrc];
            }
            else {
                arrayDst[(iRowDst*nColSigned+iColDst)*nSecSigned+iSecDst] = 0;
            }
        }
    }
#else
    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    nColSigned = (ptrdiff_t) nCol;
    ptrdiff_t    nSecSigned = (ptrdiff_t) nSec;
    ptrdiff_t    iRowDst = blockDim.z * blockIdx.z + threadIdx.z;
    ptrdiff_t    iColDst = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iSecDst = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iRowSrc, iColSrc, iSecSrc;

    if (iRowDst < nRowSigned && iColDst < nColSigned & iSecDst < nSecSigned) {
        iRowSrc = iRowDst - nRowOff;
        iColSrc = iColDst - nColOff;
        iSecSrc = iSecDst - nSecOff;

        if ((iRowSrc >= 0) && (iRowSrc < nRowSigned) &&
            (iColSrc >= 0) && (iColSrc < nColSigned) &&
            (iSecSrc >= 0) && (iSecSrc < nSecSigned)) {
            arrayDst[(iRowDst*nColSigned+iColDst)*nSecSigned+iSecDst] =
            arraySrc[(iRowSrc*nColSigned+iColSrc)*nSecSigned+iSecSrc];
        }
        else {
            arrayDst[(iRowDst*nColSigned+iColDst)*nSecSigned+iSecDst] = 0;
        }
    }
#endif
}

template <typename T>
void cuda_array_shift(const T* const arraySrc, T* const arrayDst,
                      size_t    nRow,    size_t    nCol,    size_t    nSec,
                      ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

#ifdef __GEM_CUDA_ARCH_HOST_130__
    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nSec, dimBlock.x),
                     iDivUp(nCol, dimBlock.y));
#else
    dim3    dimBlock(BLOCK_3D_NSEC, BLOCK_3D_NCOL, BLOCK_3D_NROW);
    dim3    dimGrid (iDivUp(nSec, dimBlock.x),
                     iDivUp(nCol, dimBlock.y),
                     iDivUp(nRow, dimBlock.z));
#endif

    dev_array_shift<<<dimGrid, dimBlock>>>
        (arraySrc, arrayDst, nRow, nCol, nSec, nRowOff, nColOff, nSecOff);
}

// instantiation
template
void cuda_array_shift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                                size_t    nRow,    size_t    nCol,    size_t    nSec,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_shift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                                size_t    nRow,    size_t    nCol,    size_t    nSec,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_shift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                                size_t    nRow,    size_t    nCol,    size_t    nSec,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_shift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                                size_t    nRow,    size_t    nCol,    size_t    nSec,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_shift<float   >(const float*    const arraySrc, float*    const arrayDst,
                                size_t    nRow,    size_t    nCol,    size_t    nSec,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_shift<double  >(const double*   const arraySrc, double*   const arrayDst,
                                size_t    nRow,    size_t    nCol,    size_t    nSec,
                                ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

} // namespace gem

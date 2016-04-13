/***********************************************************************
 *  File:       array_circshifting.cu
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
 * Circular shifting
 ****************************************/

// 1D
template <typename T> __global__
void dev_array_circshift(const T* const arraySrc, T* const arrayDst,
                         size_t    nRow,
                         ptrdiff_t nRowOff)
{
    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    iRowSrc = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iRowDst;

    if (iRowSrc < nRowSigned) {
        iRowDst = (iRowSrc+nRowOff) % nRowSigned;
        iRowDst = (iRowDst < 0) ? (iRowDst+nRowSigned) : (iRowDst);

        arrayDst[iRowDst] = arraySrc[iRowSrc];
    }
}

template <typename T>
void cuda_array_circshift(const T* const arraySrc, T* const arrayDst,
                          size_t    nRow,
                          ptrdiff_t nRowOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    dev_array_circshift<<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arraySrc, arrayDst, nRow, nRowOff);
}

// instantiation
template
void cuda_array_circshift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                                    size_t    nRow,
                                    ptrdiff_t nRowOff);
template
void cuda_array_circshift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                                    size_t    nRow,
                                    ptrdiff_t nRowOff);
template
void cuda_array_circshift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                                    size_t    nRow,
                                    ptrdiff_t nRowOff);
template
void cuda_array_circshift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                                    size_t    nRow,
                                    ptrdiff_t nRowOff);
template
void cuda_array_circshift<float   >(const float*    const arraySrc, float*    const arrayDst,
                                    size_t    nRow,
                                    ptrdiff_t nRowOff);
template
void cuda_array_circshift<double  >(const double*   const arraySrc, double*   const arrayDst,
                                    size_t    nRow,
                                    ptrdiff_t nRowOff);

// 2D
template <typename T> __global__
void dev_array_circshift(const T* const arraySrc, T* const arrayDst,
                         size_t    nRow,    size_t    nCol,
                         ptrdiff_t nRowOff, ptrdiff_t nColOff)
{
    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    nColSigned = (ptrdiff_t) nCol;
    ptrdiff_t    iRowSrc = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iColSrc = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iRowDst, iColDst;

    if (iRowSrc < nRowSigned && iColSrc < nColSigned) {
        iRowDst = (iRowSrc+nRowOff) % nRowSigned;
        iRowDst = (iRowDst < 0) ? (iRowDst+nRowSigned) : (iRowDst);
        iColDst = (iColSrc+nColOff) % nColSigned;
        iColDst = (iColDst < 0) ? (iColDst+nColSigned) : (iColDst);

        arrayDst[iRowDst*nColSigned+iColDst] =
        arraySrc[iRowSrc*nColSigned+iColSrc];
    }
}

template <typename T>
void cuda_array_circshift(const T* const arraySrc, T* const arrayDst,
                          size_t    nRow,    size_t    nCol,
                          ptrdiff_t nRowOff, ptrdiff_t nColOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nCol, dimBlock.x),
                     iDivUp(nRow, dimBlock.y));

    dev_array_circshift<<<dimGrid, dimBlock>>>
        (arraySrc, arrayDst, nRow, nCol, nRowOff, nColOff);
}

// instantiation
template
void cuda_array_circshift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                                    size_t    nRow,    size_t    nCol,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_circshift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                                    size_t    nRow,    size_t    nCol,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_circshift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                                    size_t    nRow,    size_t    nCol,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_circshift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                                    size_t    nRow,    size_t    nCol,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_circshift<float   >(const float*    const arraySrc, float*    const arrayDst,
                                    size_t    nRow,    size_t    nCol,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_circshift<double  >(const double*   const arraySrc, double*   const arrayDst,
                                    size_t    nRow,    size_t    nCol,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D
template <typename T> __global__
void dev_array_circshift(const T* const arraySrc, T* const arrayDst,
                         size_t    nRow,    size_t    nCol,    size_t    nSec,
                         ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    nColSigned = (ptrdiff_t) nCol;
    ptrdiff_t    nSecSigned = (ptrdiff_t) nSec;
    ptrdiff_t    iRowSrc;
    ptrdiff_t    iColSrc = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iSecSrc = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iRowDst, iColDst, iSecDst;

    if (iColSrc < nColSigned & iSecSrc < nSecSigned) {
        iColDst = (iColSrc+nColOff) % nColSigned;
        iColDst = (iColDst < 0) ? (iColDst+nColSigned) : (iColDst);
        iSecDst = (iSecSrc+nSecOff) % nSecSigned;
        iSecDst = (iSecDst < 0) ? (iSecDst+nSecSigned) : (iSecDst);

        for (iRowSrc = 0; iRowSrc < nRowSigned; iRowSrc++) {
            iRowDst = (iRowSrc+nRowOff) % nRowSigned;
            iRowDst = (iRowDst < 0) ? (iRowDst+nRowSigned) : (iRowDst);

            arrayDst[(iRowDst*nColSigned+iColDst)*nSecSigned+iSecDst] =
            arraySrc[(iRowSrc*nColSigned+iColSrc)*nSecSigned+iSecSrc];
        }
    }
#else
    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    nColSigned = (ptrdiff_t) nCol;
    ptrdiff_t    nSecSigned = (ptrdiff_t) nSec;
    ptrdiff_t    iRowSrc = blockDim.z * blockIdx.z + threadIdx.z;
    ptrdiff_t    iColSrc = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iSecSrc = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iRowDst, iColDst, iSecDst;

    if (iRowSrc < nRowSigned && iColSrc < nColSigned & iSecSrc < nSecSigned) {
        iRowDst = (iRowSrc+nRowOff) % nRowSigned;
        iRowDst = (iRowDst < 0) ? (iRowDst+nRowSigned) : (iRowDst);
        iColDst = (iColSrc+nColOff) % nColSigned;
        iColDst = (iColDst < 0) ? (iColDst+nColSigned) : (iColDst);
        iSecDst = (iSecSrc+nSecOff) % nSecSigned;
        iSecDst = (iSecDst < 0) ? (iSecDst+nSecSigned) : (iSecDst);

        arrayDst[(iRowDst*nColSigned+iColDst)*nSecSigned+iSecDst] =
        arraySrc[(iRowSrc*nColSigned+iColSrc)*nSecSigned+iSecSrc];
    }
#endif
}

template <typename T>
void cuda_array_circshift(const T* const arraySrc, T* const arrayDst,
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

    dev_array_circshift<<<dimGrid, dimBlock>>>
        (arraySrc, arrayDst, nRow, nCol, nSec, nRowOff, nColOff, nSecOff);
}

// instantiation
template
void cuda_array_circshift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                                    size_t    nRow,    size_t    nCol,    size_t    nSec,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_circshift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                                    size_t    nRow,    size_t    nCol,    size_t    nSec,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_circshift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                                    size_t    nRow,    size_t    nCol,    size_t    nSec,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_circshift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                                    size_t    nRow,    size_t    nCol,    size_t    nSec,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_circshift<float   >(const float*    const arraySrc, float*    const arrayDst,
                                    size_t    nRow,    size_t    nCol,    size_t    nSec,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_circshift<double  >(const double*   const arraySrc, double*   const arrayDst,
                                    size_t    nRow,    size_t    nCol,    size_t    nSec,
                                    ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

} // namespace gem

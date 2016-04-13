/***********************************************************************
 *  File:       array_circshifting.cpp
 *
 *  Purpose:    Implementation of array-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "array.hpp"

namespace gem {

/*****************************************
 * Circular shifting
 ****************************************/

// 1D
template <typename T>
void array_circshift(const T* const arraySrc, T* const arrayDst,
                     size_t    nRow,
                     ptrdiff_t nRowOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    iRowSrc;
    ptrdiff_t    iRowDst;

    #pragma omp parallel for private(iRowSrc,iRowDst)
    for (iRowSrc = 0; iRowSrc < nRowSigned; iRowSrc++) {
        iRowDst = (iRowSrc+nRowOff) % nRowSigned;
        iRowDst = (iRowDst < 0) ? (iRowDst+nRowSigned) : (iRowDst);

        arrayDst[iRowDst] = arraySrc[iRowSrc];
    }
}

// instantiation
template
void array_circshift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                               size_t    nRow,
                               ptrdiff_t nRowOff);
template
void array_circshift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                               size_t    nRow,
                               ptrdiff_t nRowOff);
template
void array_circshift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                               size_t    nRow,
                               ptrdiff_t nRowOff);
template
void array_circshift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                               size_t    nRow,
                               ptrdiff_t nRowOff);
template
void array_circshift<float   >(const float*    const arraySrc, float*    const arrayDst,
                               size_t    nRow,
                               ptrdiff_t nRowOff);
template
void array_circshift<double  >(const double*   const arraySrc, double*   const arrayDst,
                               size_t    nRow,
                               ptrdiff_t nRowOff);

// 2D
template <typename T>
void array_circshift(const T* const arraySrc, T* const arrayDst,
                     size_t    nRow,    size_t    nCol,
                     ptrdiff_t nRowOff, ptrdiff_t nColOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);

    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    nColSigned = (ptrdiff_t) nCol;
    ptrdiff_t    iSrc1, iRowSrc, iColSrc;
    ptrdiff_t    iDst1, iRowDst, iColDst;

    #pragma omp parallel for private(iRowSrc,iRowDst,iSrc1,iDst1,iColSrc,iColDst)
    for (iRowSrc = 0; iRowSrc < nRowSigned; iRowSrc++) {
        iRowDst = (iRowSrc+nRowOff) % nRowSigned;
        iRowDst = (iRowDst < 0) ? (iRowDst+nRowSigned) : (iRowDst);
        iSrc1 = iRowSrc*nColSigned;
        iDst1 = iRowDst*nColSigned;

        for (iColSrc = 0; iColSrc < nColSigned; iColSrc++) {
            iColDst = (iColSrc+nColOff) % nColSigned;
            iColDst = (iColDst < 0) ? (iColDst+nColSigned) : (iColDst);

            arrayDst[iDst1+iColDst] = arraySrc[iSrc1+iColSrc];
        }
    }
}

// instantiation
template
void array_circshift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                               size_t    nRow,    size_t    nCol,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_circshift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                               size_t    nRow,    size_t    nCol,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_circshift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                               size_t    nRow,    size_t    nCol,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_circshift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                               size_t    nRow,    size_t    nCol,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_circshift<float   >(const float*    const arraySrc, float*    const arrayDst,
                             size_t    nRow,    size_t    nCol,
                             ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_circshift<double  >(const double*   const arraySrc, double*   const arrayDst,
                               size_t    nRow,    size_t    nCol,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D
template <typename T>
void array_circshift(const T* const arraySrc, T* const arrayDst,
                     size_t    nRow,    size_t    nCol,    size_t    nSec,
                     ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    nColSigned = (ptrdiff_t) nCol;
    ptrdiff_t    nSecSigned = (ptrdiff_t) nSec;
    ptrdiff_t    iSrc1, iSrc2, iRowSrc, iColSrc, iSecSrc;
    ptrdiff_t    iDst1, iDst2, iRowDst, iColDst, iSecDst;
    ptrdiff_t    nColSecSigned = (ptrdiff_t) (nCol * nSec);

    #pragma omp parallel for private(iRowSrc,iRowDst,iSrc1,iDst1,iColSrc,iColDst,iSrc2,iDst2,iSecSrc,iSecDst)
    for (iRowSrc = 0; iRowSrc < nRowSigned; iRowSrc++) {
        iRowDst = (iRowSrc+nRowOff) % nRowSigned;
        iRowDst = (iRowDst < 0) ? (iRowDst+nRowSigned) : (iRowDst);
        iSrc1 = iRowSrc*nColSecSigned;
        iDst1 = iRowDst*nColSecSigned;

        for (iColSrc = 0; iColSrc < nColSigned; iColSrc++) {
            iColDst = (iColSrc+nColOff) % nColSigned;
            iColDst = (iColDst < 0) ? (iColDst+nColSigned) : (iColDst);
            iSrc2 = iSrc1 + iColSrc*nSecSigned;
            iDst2 = iDst1 + iColDst*nSecSigned;

            for (iSecSrc = 0; iSecSrc < nSecSigned; iSecSrc++) {
                iSecDst = (iSecSrc+nSecOff) % nSecSigned;
                iSecDst = (iSecDst < 0) ? (iSecDst+nSecSigned) : (iSecDst);

                arrayDst[iDst2+iSecDst] = arraySrc[iSrc2+iSecSrc];
            }
        }
    }
}

// instantiation
template
void array_circshift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                               size_t    nRow,    size_t    nCol,    size_t    nSec,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_circshift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                               size_t    nRow,    size_t    nCol,    size_t    nSec,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_circshift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                               size_t    nRow,    size_t    nCol,    size_t    nSec,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_circshift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                               size_t    nRow,    size_t    nCol,    size_t    nSec,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_circshift<float   >(const float*    const arraySrc, float*    const arrayDst,
                               size_t    nRow,    size_t    nCol,    size_t    nSec,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_circshift<double  >(const double*   const arraySrc, double*   const arrayDst,
                               size_t    nRow,    size_t    nCol,    size_t    nSec,
                               ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

} // namespace gem

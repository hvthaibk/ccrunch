/***********************************************************************
 *  File:       array_shifting.cpp
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
 * Shifting
 ****************************************/

// 1D
template <typename T>
void array_shift(const T* const arraySrc, T* const arrayDst,
                 size_t    nRow,
                 ptrdiff_t nRowOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    ptrdiff_t    nRowSigned = (ptrdiff_t) nRow;
    ptrdiff_t    iRowSrc;
    ptrdiff_t    iRowDst;

    #pragma omp parallel for private(iRowDst,iRowSrc)
    for (iRowDst = 0; iRowDst < nRowSigned; iRowDst++) {
        iRowSrc = iRowDst - nRowOff;

        if ((iRowSrc >= 0) && (iRowSrc < nRowSigned)) {
            arrayDst[iRowDst] = arraySrc[iRowSrc];
        }
        else {
            arrayDst[iRowDst] = 0;
        }
    }
}

// instantiation
template
void array_shift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                           size_t    nRow,
                           ptrdiff_t nRowOff);
template
void array_shift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                           size_t    nRow,
                           ptrdiff_t nRowOff);
template
void array_shift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                           size_t    nRow,
                           ptrdiff_t nRowOff);
template
void array_shift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                           size_t    nRow,
                           ptrdiff_t nRowOff);
template
void array_shift<float   >(const float*    const arraySrc, float*    const arrayDst,
                           size_t    nRow,
                           ptrdiff_t nRowOff);
template
void array_shift<double  >(const double*   const arraySrc, double*   const arrayDst,
                           size_t    nRow,
                           ptrdiff_t nRowOff);

// 2D
template <typename T>
void array_shift(const T* const arraySrc, T* const arrayDst,
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

    #pragma omp parallel for private(iRowDst,iRowSrc,iSrc1,iDst1,iColDst,iColSrc)
    for (iRowDst = 0; iRowDst < nRowSigned; iRowDst++) {
        iRowSrc = iRowDst - nRowOff;
        iSrc1 = iRowSrc*nColSigned;
        iDst1 = iRowDst*nColSigned;

        for (iColDst = 0; iColDst < nColSigned; iColDst++) {
            iColSrc = iColDst - nColOff;

            if ((iRowSrc >= 0) && (iRowSrc < nRowSigned) &&
                (iColSrc >= 0) && (iColSrc < nColSigned)) {
                arrayDst[iDst1+iColDst] = arraySrc[iSrc1+iColSrc];
            }
            else {
                arrayDst[iDst1+iColDst] = 0;
            }
        }
    }
}

// instantiation
template
void array_shift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                           size_t    nRow,    size_t    nCol,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_shift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                           size_t    nRow,    size_t    nCol,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_shift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                           size_t    nRow,    size_t    nCol,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_shift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                           size_t    nRow,    size_t    nCol,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_shift<float   >(const float*    const arraySrc, float*    const arrayDst,
                           size_t    nRow,    size_t    nCol,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_shift<double  >(const double*   const arraySrc, double*   const arrayDst,
                           size_t    nRow,    size_t    nCol,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D
template <typename T>
void array_shift(const T* const arraySrc, T* const arrayDst,
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

    #pragma omp parallel for private(iRowDst,iRowSrc,iSrc1,iDst1,iColDst,iColSrc,iSrc2,iDst2,iSecDst,iSecSrc)
    for (iRowDst = 0; iRowDst < nRowSigned; iRowDst++) {
        iRowSrc = iRowDst - nRowOff;
        iSrc1 = iRowSrc*nColSecSigned;
        iDst1 = iRowDst*nColSecSigned;

        for (iColDst = 0; iColDst < nColSigned; iColDst++) {
            iColSrc = iColDst - nColOff;
            iSrc2 = iSrc1 + iColSrc*nSecSigned;
            iDst2 = iDst1 + iColDst*nSecSigned;

            for (iSecDst = 0; iSecDst < nSecSigned; iSecDst++) {
                iSecSrc = iSecDst - nSecOff;

                if ((iRowSrc >= 0) && (iRowSrc < nRowSigned) &&
                    (iColSrc >= 0) && (iColSrc < nColSigned) &&
                    (iSecSrc >= 0) && (iSecSrc < nSecSigned)) {

                    arrayDst[iDst2+iSecDst] = arraySrc[iSrc2+iSecSrc];
                }
                else {
                    arrayDst[iDst2+iSecDst] = 0;
                }
            }
        }
    }
}

// instantiation
template
void array_shift<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                           size_t    nRow,    size_t    nCol,    size_t    nSec,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_shift<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                           size_t    nRow,    size_t    nCol,    size_t    nSec,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_shift<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                           size_t    nRow,    size_t    nCol,    size_t    nSec,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_shift<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                           size_t    nRow,    size_t    nCol,    size_t    nSec,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_shift<float   >(const float*    const arraySrc, float*    const arrayDst,
                           size_t    nRow,    size_t    nCol,    size_t    nSec,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_shift<double  >(const double*   const arraySrc, double*   const arrayDst,
                           size_t    nRow,    size_t    nCol,    size_t    nSec,
                           ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

} // namespace gem

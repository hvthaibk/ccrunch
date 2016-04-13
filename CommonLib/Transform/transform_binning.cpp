/***********************************************************************
 *  File:       transform_binning.cpp
 *
 *  Purpose:    Implementation of transformation functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "transform.hpp"

namespace gem {

/**************************
 * Binning
 *************************/

// 1D
template <typename T>
void transform_bin(const T* const arraySrc, T* const arrayDst,
                   size_t nRow,
                   size_t factor)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0);
    assert(factor > 0);

    size_t    iRowSrc;
    size_t    nRowDst;
    T         sumTmp;

    nRowDst = transform_binSize(nRow, factor);

    for (size_t iRowDst = 0; iRowDst < nRowDst; iRowDst++) {
        iRowSrc = factor * iRowDst;

        sumTmp = 0;

        for (size_t iRow = iRowSrc; iRow < iRowSrc+factor; iRow++) {
            sumTmp += arraySrc[iRow];
        }

        arrayDst[iRowDst] = sumTmp / (T) factor;
    }
}

// instantiation
template
void transform_bin<float >(const float*  const arraySrc, float*  const arrayDst,
                           size_t nRow,
                           size_t factor);
template
void transform_bin<double>(const double* const arraySrc, double* const arrayDst,
                           size_t nRow,
                           size_t factor);

// 2D
template <typename T>
void transform_bin(const T* const arraySrc, T* const arrayDst,
                   size_t nRow, size_t nCol,
                   size_t factor)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);
    assert(factor > 0);

    size_t    iRowSrc, iColSrc;
    size_t    nRowDst, nColDst;
    T         sumTmp;

    nRowDst = transform_binSize(nRow, factor);
    nColDst = transform_binSize(nCol, factor);

    for (size_t iRowDst = 0; iRowDst < nRowDst; iRowDst++) {
        iRowSrc = factor * iRowDst;

        for (size_t iColDst = 0; iColDst < nColDst; iColDst++) {
            iColSrc = factor * iColDst;

            sumTmp = 0;

            for (size_t iRow = iRowSrc; iRow < iRowSrc+factor; iRow++) {
                for (size_t iCol = iColSrc; iCol < iColSrc+factor; iCol++) {
                    sumTmp += arraySrc[iRow*nCol+iCol];
                }
            }

            arrayDst[iRowDst*nColDst+iColDst] = sumTmp / (T) pow2(factor);
        }
    }
}

// instantiation
template
void transform_bin<float >(const float*  const arraySrc, float*  const arrayDst,
                           size_t nRow, size_t nCol,
                           size_t factor);
template
void transform_bin<double>(const double* const arraySrc, double* const arrayDst,
                           size_t nRow, size_t nCol,
                           size_t factor);

// 3D
template <typename T>
void transform_bin(const T* const arraySrc, T* const arrayDst,
                   size_t nRow, size_t nCol, size_t nSec,
                   size_t factor)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);
    assert(factor > 0);

    size_t    iRowSrc, iColSrc, iSecSrc;
    size_t    nRowDst, nColDst, nSecDst;
    T         sumTmp;

    nRowDst = transform_binSize(nRow, factor);
    nColDst = transform_binSize(nCol, factor);
    nSecDst = transform_binSize(nSec, factor);

    for (size_t iRowDst = 0; iRowDst < nRowDst; iRowDst++) {
        iRowSrc = factor * iRowDst;

        for (size_t iColDst = 0; iColDst < nColDst; iColDst++) {
            iColSrc = factor * iColDst;

            for (size_t iSecDst = 0; iSecDst < nSecDst; iSecDst++) {
                iSecSrc = factor*iSecDst;

                sumTmp = 0;

                for (size_t iRow = iRowSrc; iRow < iRowSrc+factor; iRow++) {
                    for (size_t iCol = iColSrc; iCol < iColSrc+factor; iCol++) {
                        for (size_t iSec = iSecSrc; iSec < iSecSrc+factor; iSec++) {
                            sumTmp += arraySrc[iRow*nCol*nSec+iCol*nSec+iSec];
                        }
                    }
                }

                arrayDst[iRowDst*nColDst*nSecDst+iColDst*nSecDst+iSecDst] = sumTmp / (T) pow3(factor);
            }
        }
    }
}

// instantiation
template
void transform_bin<float >(const float*  const arraySrc, float*  const arrayDst,
                           size_t nRow, size_t nCol, size_t nSec,
                           size_t factor);
template
void transform_bin<double>(const double* const arraySrc, double* const arrayDst,
                           size_t nRow, size_t nCol, size_t nSec,
                           size_t factor);

} // namespace gem

/***********************************************************************
 *  File:       transform_dilation.cpp
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

#include "filter.hpp"
#include "xcorr.hpp"

namespace gem {

/**************************
 * Dilation
 *************************/

// 1D
template <typename T>
void transform_dilate(const T* const arraySrc, size_t nRow,
                            T* const arrayDst,
                            T        radius)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 2);
    assert(radius > 0);

    // create structuring element
    T*          seData = NULL;
    size_t      seRad  = (size_t) std::ceil(radius);
    size_t      nRowSE = 2 * seRad + 1;

    array_new  (seData, nRowSE);
    filter_disk(seData, nRowSE, radius, 0);

    // extract boundary points
    T           bpFilter[3] = {-1, 2, -1};
    T*          bpData = NULL;
    size_t      nRowBP = nRow-2;

    array_new(bpData, nRowBP);

    xcorr(bpFilter, 3,
          arraySrc, nRow,
          bpData,
          XCORR_RES_VALID);

    // scan boundary points with structuring element
    array_memcpy(arrayDst, arraySrc, nRow);

    for (size_t iRow = 1, idxBP = 0; iRow < nRow-1; iRow++, idxBP++) {

        if (bpData[idxBP] > (T) 0.5) {

            for (size_t iRowSE = 0, idxSE = 0; iRowSE < nRowSE; iRowSE++, idxSE++) {

                if (seData[idxSE] > 0 &&
                    iRow+iRowSE >= seRad && iRow+iRowSE < nRow+seRad) {

                    arrayDst[iRow+iRowSE-seRad] = 1;
                }
            }
        }
    }

    // clean memory
    array_delete(seData);
    array_delete(bpData);
}

// instantiation
template
void transform_dilate<float >(const float*  const arraySrc, size_t nRow,
                                    float*  const arrayDst,
                                    float         radius);
template
void transform_dilate<double>(const double* const arraySrc, size_t nRow,
                                    double* const arrayDst,
                                    double        radius);

// 2D
template <typename T>
void transform_dilate(const T* const arraySrc, size_t nRow, size_t nCol,
                            T* const arrayDst,
                            T        radius)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 2 && nCol > 2);
    assert(radius > 0);

    // create structuring element
    T*          seData = NULL;
    size_t      seRad  = (size_t) std::ceil(radius);
    size_t      nRowSE = 2 * seRad + 1;
    size_t      nColSE = 2 * seRad + 1;

    array_new  (seData, nRowSE*nColSE);
    filter_disk(seData, nRowSE, nColSE, radius, 0);

    // extract boundary points
    T           bpFilter[9] = {-1, -1, -1, -1, 8, -1, -1, -1, -1};
    T*          bpData = NULL;
    size_t      nRowBP = nRow-2;
    size_t      nColBP = nCol-2;

    array_new(bpData, nRowBP*nColBP);

    xcorr(bpFilter, 3, 3,
          arraySrc, nRow, nCol,
          bpData,
          XCORR_RES_VALID);



    // scan boundary points with structuring element
    array_memcpy(arrayDst, arraySrc, nRow*nCol);

    for (size_t iRow = 1, idxBP = 0; iRow < nRow-1; iRow++) {
        for (size_t iCol = 1; iCol < nCol-1; iCol++, idxBP++) {

            if (bpData[idxBP] > (T) 0.5) {

                for (size_t iRowSE = 0, idxSE = 0; iRowSE < nRowSE; iRowSE++) {
                    for (size_t iColSE = 0; iColSE < nColSE; iColSE++, idxSE++) {

                        if (seData[idxSE] > 0 &&
                            iRow+iRowSE >= seRad && iRow+iRowSE < nRow+seRad &&
                            iCol+iColSE >= seRad && iCol+iColSE < nCol+seRad) {

                                arrayDst[(iRow+iRowSE-seRad)*nCol+(iCol+iColSE-seRad)] = 1;
                        }
                    }
                }
            }
        }
    }

    // clean memory
    array_delete(seData);
    array_delete(bpData);
}

// instantiation
template
void transform_dilate<float >(const float*  const arraySrc, size_t nRow, size_t nCol,
                                    float*  const arrayDst,
                                    float         radius);
template
void transform_dilate<double>(const double* const arraySrc, size_t nRow, size_t nCol,
                                    double* const arrayDst,
                                    double        radius);

// 3D
template <typename T>
void transform_dilate(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                            T* const arrayDst,
                            T        radius)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 2 && nCol > 2 && nSec > 2);
    assert(radius > 0);

    // create structuring element
    T*          seData = NULL;
    size_t      seRad  = (size_t) std::ceil(radius);
    size_t      nRowSE = 2 * seRad + 1;
    size_t      nColSE = 2 * seRad + 1;
    size_t      nSecSE = 2 * seRad + 1;

    array_new  (seData, nRowSE*nColSE*nSecSE);
    filter_disk(seData, nRowSE, nColSE, nSecSE, radius, 0);

    // extract boundary points
    T           bpFilter[27] = {-1, -1, -1, -1, -1, -1, -1, -1, -1,
                                -1, -1, -1, -1, 26, -1, -1, -1, -1,
                                -1, -1, -1, -1, -1, -1, -1, -1, -1};
    T*          bpData = NULL;
    size_t      nRowBP = nRow-2;
    size_t      nColBP = nCol-2;
    size_t      nSecBP = nSec-2;

    array_new(bpData, nRowBP*nColBP*nSecBP);

    xcorr(bpFilter, 3, 3, 3,
          arraySrc, nRow, nCol, nSec,
          bpData,
          XCORR_RES_VALID);

    // scan boundary points with structuring element
    array_memcpy(arrayDst, arraySrc, nRow*nCol*nSec);

    for (size_t iRow = 1, idxBP = 0; iRow < nRow-1; iRow++) {
        for (size_t iCol = 1; iCol < nCol-1; iCol++) {
            for (size_t iSec = 1; iSec < nSec-1; iSec++, idxBP++) {

                if (bpData[idxBP] > (T) 0.5) {

                    for (size_t iRowSE = 0, idxSE = 0; iRowSE < nRowSE; iRowSE++) {
                        for (size_t iColSE = 0; iColSE < nColSE; iColSE++) {
                            for (size_t iSecSE = 0; iSecSE < nSecSE; iSecSE++, idxSE++) {

                                if (seData[idxSE] > 0 &&
                                    iRow+iRowSE >= seRad && iRow+iRowSE < nRow+seRad &&
                                    iCol+iColSE >= seRad && iCol+iColSE < nCol+seRad &&
                                    iSec+iSecSE >= seRad && iSec+iSecSE < nSec+seRad) {

                                        arrayDst[(iRow+iRowSE-seRad)*nCol*nSec+(iCol+iColSE-seRad)*nSec+(iSec+iSecSE-seRad)] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // clean memory
    array_delete(seData);
    array_delete(bpData);
}

// instantiation
template
void transform_dilate<float >(const float*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                    float*  const arrayDst,
                                    float         radius);
template
void transform_dilate<double>(const double* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                    double* const arrayDst,
                                    double        radius);


} // namespace gem

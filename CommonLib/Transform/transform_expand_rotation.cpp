/***********************************************************************
 *  File:       transform_expand_rotation.cpp
 *
 *  Purpose:    Implementation of expansion functions
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
 * Expansion
 *************************/

// 2D
template <typename T>
void transform_expand_rotate_valid(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                         T*&      arrayDst, size_t& nRowDst, size_t& nColDst)
{
    assert(arraySrc != NULL && arrayDst == NULL);
    assert(nRowSrc > 0 && nColSrc > 0);

    T         cRow, cCol, cOff = 0.5;
    T         pixRow, pixRow1, pixRow2, pixRowMax,
              pixCol, pixCol1, pixCol2, pixColMax;
    T         dist[4], distToCenter, distToCenterMax = 0;
    size_t    diagLength;

    T         *arrayCrp = NULL;
    size_t    nRowCrp, nColCrp;
    size_t    nRowSrcOff, nColSrcOff;
    size_t    nRowDstOff, nColDstOff;

    // centroid
    pixRowMax = cRow = (T) (nRowSrc-1) / 2;
    pixColMax = cCol = (T) (nColSrc-1) / 2;

    for (size_t i = 0; i < nRowSrc*nColSrc; i++) {

        if (arraySrc[i] > GEM_FLOATING_TOLERANCE) {

            // pixel coordinates
            pixRow = std::floor((T) i / (T) nColSrc);
            pixCol = (T) i - pixRow * (T) nColSrc;

            if (std::abs(pixRow-cRow) > std::abs(pixRowMax-cRow)) {
                pixRowMax = pixRow;
            }
            if (std::abs(pixCol-cCol) > std::abs(pixColMax-cCol)) {
                pixColMax = pixCol;
            }

            pixRow1 = std::abs(pixRow-cRow-cOff);
            pixRow2 = std::abs(pixRow-cRow+cOff);
            pixCol1 = std::abs(pixCol-cCol-cOff);
            pixCol2 = std::abs(pixCol-cCol+cOff);

            dist[0] = std::sqrt(pow2(pixRow1)+pow2(pixCol1));
            dist[1] = std::sqrt(pow2(pixRow1)+pow2(pixCol2));
            dist[2] = std::sqrt(pow2(pixRow2)+pow2(pixCol1));
            dist[3] = std::sqrt(pow2(pixRow2)+pow2(pixCol2));
            distToCenter = array_reduce_max(dist, 4);

            if (distToCenterMax < distToCenter) {
                distToCenterMax = distToCenter;
            }
        }
    }

    diagLength = (size_t) std::ceil(2*distToCenterMax);
    require(diagLength > 0, "transform_expand_rotate_valid(): input array contains only zero values");

    // crop array
    nRowSrcOff = (size_t) ((pixRowMax < cRow) ? pixRowMax : ((T) nRowSrc-1-pixRowMax));
    nColSrcOff = (size_t) ((pixColMax < cCol) ? pixColMax : ((T) nColSrc-1-pixColMax));
    nRowCrp = nRowSrc - 2*nRowSrcOff;
    nColCrp = nColSrc - 2*nColSrcOff;

    array_new(arrayCrp, nRowCrp*nColCrp);
    array_crop(arraySrc, nRowSrc,    nColSrc,
               arrayCrp, nRowCrp,    nColCrp,
                         nRowSrcOff, nColSrcOff);

    // dest array
    nRowDst = ((diagLength-nRowSrc) % 2) ? diagLength+1 : diagLength;
    nColDst = ((diagLength-nColSrc) % 2) ? diagLength+1 : diagLength;
    nRowDst = std::max(nRowDst, nRowCrp);
    nColDst = std::max(nColDst, nColCrp);
    nRowDstOff = (nRowDst-nRowCrp) / 2;
    nColDstOff = (nColDst-nColCrp) / 2;

    array_new_zero(arrayDst, nRowDst*nColDst);

    array_replace<T,T>(arrayDst, nRowDst,    nColDst,
                       arrayCrp, nRowCrp,    nColCrp,
                                 nRowDstOff, nColDstOff);

    array_delete(arrayCrp);
}

template <typename T>
void transform_expand_rotate_full(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                        T*&      arrayDst, size_t& nRowDst, size_t& nColDst)
{
    assert(arraySrc != NULL && arrayDst == NULL);
    assert(nRowSrc > 0 && nColSrc > 0);

    size_t    diagLength;
    size_t    nRowOff, nColOff;

    diagLength = (size_t) std::ceil(dist((T)0,(T)0,(T)nRowSrc,(T)nColSrc));

    nRowDst = ((diagLength-nRowSrc) % 2) ? diagLength+1 : diagLength;
    nColDst = ((diagLength-nColSrc) % 2) ? diagLength+1 : diagLength;
    //nRowDst = diagLength;
    //nColDst = diagLength;

    nRowOff = (nRowDst-nRowSrc) / 2;
    nColOff = (nColDst-nColSrc) / 2;

    array_new_zero(arrayDst, nRowDst*nColDst);

    array_replace<T,T>(arrayDst, nRowDst, nColDst,
                       arraySrc, nRowSrc, nColSrc,
                                 nRowOff, nColOff);
}

// instantiation
template
void transform_expand_rotate_valid<float >(const float*  const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                                 float*&       arrayDst, size_t& nRowDst, size_t& nColDst);
template
void transform_expand_rotate_valid<double>(const double* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                                 double*&      arrayDst, size_t& nRowDst, size_t& nColDst);
template
void transform_expand_rotate_full <float >(const float*  const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                                 float*&       arrayDst, size_t& nRowDst, size_t& nColDst);
template
void transform_expand_rotate_full <double>(const double* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                                 double*&      arrayDst, size_t& nRowDst, size_t& nColDst);

// 3D
template <typename T>
void transform_expand_rotate_valid(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                         T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst)
{
    assert(arraySrc != NULL && arrayDst == NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);

    T         cRow, cCol, cSec, cOff = 0.5;
    T         pixRow, pixRow1, pixRow2, pixRowMax,
              pixCol, pixCol1, pixCol2, pixColMax,
              pixSec, pixSec1, pixSec2, pixSecMax;
    T         dist[8], distToCenter, distToCenterMax = 0;
    size_t    diagLength;

    T         *arrayCrp = NULL;
    size_t    nRowCrp, nColCrp, nSecCrp;
    size_t    nRowSrcOff, nColSrcOff, nSecSrcOff;
    size_t    nRowDstOff, nColDstOff, nSecDstOff;

    // centroid
    pixRowMax = cRow = (T) (nRowSrc-1) / 2;
    pixColMax = cCol = (T) (nColSrc-1) / 2;
    pixSecMax = cSec = (T) (nSecSrc-1) / 2;

    for (size_t i = 0; i < nRowSrc*nColSrc*nSecSrc; i++) {

        if (arraySrc[i] > GEM_FLOATING_TOLERANCE) {

            // pixel coordinates
            T     iTmp = (T) i;
            pixRow = std::floor(iTmp / (T) (nColSrc*nSecSrc));
            iTmp   = (T) i - pixRow * (T) (nColSrc*nSecSrc);
            pixCol = std::floor(iTmp / (T) nSecSrc);
            pixSec = iTmp - pixCol * (T) nSecSrc;

            if (std::abs(pixRow-cRow) > std::abs(pixRowMax-cRow)) {
                pixRowMax = pixRow;
            }
            if (std::abs(pixCol-cCol) > std::abs(pixColMax-cCol)) {
                pixColMax = pixCol;
            }
            if (std::abs(pixSec-cSec) > std::abs(pixSecMax-cSec)) {
                pixSecMax = pixSec;
            }

            pixRow1 = std::abs(pixRow-cRow-cOff);
            pixRow2 = std::abs(pixRow-cRow+cOff);
            pixCol1 = std::abs(pixCol-cCol-cOff);
            pixCol2 = std::abs(pixCol-cCol+cOff);
            pixSec1 = std::abs(pixSec-cSec-cOff);
            pixSec2 = std::abs(pixSec-cSec+cOff);

            dist[0] = std::sqrt(pow2(pixRow1)+pow2(pixCol1)+pow2(pixSec1));
            dist[1] = std::sqrt(pow2(pixRow1)+pow2(pixCol1)+pow2(pixSec2));
            dist[2] = std::sqrt(pow2(pixRow1)+pow2(pixCol2)+pow2(pixSec1));
            dist[3] = std::sqrt(pow2(pixRow1)+pow2(pixCol2)+pow2(pixSec2));
            dist[4] = std::sqrt(pow2(pixRow2)+pow2(pixCol1)+pow2(pixSec1));
            dist[5] = std::sqrt(pow2(pixRow2)+pow2(pixCol1)+pow2(pixSec2));
            dist[6] = std::sqrt(pow2(pixRow2)+pow2(pixCol2)+pow2(pixSec1));
            dist[7] = std::sqrt(pow2(pixRow2)+pow2(pixCol2)+pow2(pixSec2));
            distToCenter = array_reduce_max(dist, 8);

            if (distToCenterMax < distToCenter) {
                distToCenterMax = distToCenter;
            }
        }
    }

    diagLength = (size_t) std::ceil(2*distToCenterMax);
    require(diagLength > 0, "transform_expand_rotate_valid(): input array contains only zero values");

    // crop array
    nRowSrcOff = (size_t) ((pixRowMax < cRow) ? pixRowMax : ((T) nRowSrc-1-pixRowMax));
    nColSrcOff = (size_t) ((pixColMax < cCol) ? pixColMax : ((T) nColSrc-1-pixColMax));
    nSecSrcOff = (size_t) ((pixSecMax < cSec) ? pixSecMax : ((T) nSecSrc-1-pixSecMax));
    nRowCrp = nRowSrc - 2*nRowSrcOff;
    nColCrp = nColSrc - 2*nColSrcOff;
    nSecCrp = nSecSrc - 2*nSecSrcOff;

    array_new(arrayCrp, nRowCrp*nColCrp*nSecCrp);
    array_crop(arraySrc, nRowSrc,    nColSrc,    nSecSrc,
               arrayCrp, nRowCrp,    nColCrp,    nSecCrp,
                         nRowSrcOff, nColSrcOff, nSecSrcOff);

    // dest array
    nRowDst = ((diagLength-nRowSrc) % 2) ? diagLength+1 : diagLength;
    nColDst = ((diagLength-nColSrc) % 2) ? diagLength+1 : diagLength;
    nSecDst = ((diagLength-nSecSrc) % 2) ? diagLength+1 : diagLength;
    nRowDst = std::max(nRowDst, nRowCrp);
    nColDst = std::max(nColDst, nColCrp);
    nSecDst = std::max(nSecDst, nSecCrp);
    nRowDstOff = (nRowDst-nRowCrp) / 2;
    nColDstOff = (nColDst-nColCrp) / 2;
    nSecDstOff = (nSecDst-nSecCrp) / 2;

    array_new_zero(arrayDst, nRowDst*nColDst*nSecDst);

    array_replace<T,T>(arrayDst, nRowDst,    nColDst,    nSecDst,
                       arrayCrp, nRowCrp,    nColCrp,    nSecCrp,
                                 nRowDstOff, nColDstOff, nSecDstOff);

    array_delete(arrayCrp);
}

template <typename T>
void transform_expand_rotate_full(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                        T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst)
{
    assert(arraySrc != NULL && arrayDst == NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);

    size_t    diagLength;
    size_t    nRowOff, nColOff, nSecOff;

    diagLength = (size_t) std::ceil(dist((T)0,(T)0,(T)0,(T)nRowSrc,(T)nColSrc,(T)nSecSrc));

    nRowDst = ((diagLength-nRowSrc) % 2) ? diagLength+1 : diagLength;
    nColDst = ((diagLength-nColSrc) % 2) ? diagLength+1 : diagLength;
    nSecDst = ((diagLength-nSecSrc) % 2) ? diagLength+1 : diagLength;
    //nRowDst = diagLength;
    //nColDst = diagLength;
    //nSecDst = diagLength;

    nRowOff = (nRowDst-nRowSrc) / 2;
    nColOff = (nColDst-nColSrc) / 2;
    nSecOff = (nSecDst-nSecSrc) / 2;

    array_new_zero(arrayDst, nRowDst*nColDst*nSecDst);

    array_replace<T,T>(arrayDst, nRowDst, nColDst, nSecDst,
                       arraySrc, nRowSrc, nColSrc, nSecSrc,
                                 nRowOff, nColOff, nSecOff);
}

// instantiation
template
void transform_expand_rotate_valid<float >(const float*  const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                                 float*&       arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst);
template
void transform_expand_rotate_valid<double>(const double* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                                 double*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst);
template
void transform_expand_rotate_full <float >(const float*  const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                                 float*&       arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst);
template
void transform_expand_rotate_full <double>(const double* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                                 double*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst);

} // namespace gem

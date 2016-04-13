/***********************************************************************
 *  File:       transform_translation.cpp
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
 * Translation
 *************************/

// 1D
template <typename T>
void transform_translate(const T* const arraySrc, T* const arrayDst,
                         size_t nRow,
                         T      nRowOff,
                         eInter inter)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0);

    ptrdiff_t    intRow, intRow1;
    bool         bRow, bRow1;
    T            pixRow, mapRow;
    T            v0, v1;
    T            fracRow;

    //
    // do not use OpenMP here, it slows down the program
    //
    for (size_t i = 0; i < nRow; i++) {
        // compute original pixel's coordinates
        pixRow = (T) i;

        // compute mapped pixel's coordinates
        mapRow = pixRow - nRowOff;

        switch (inter) {
            case INTER_NEAREST:
                intRow = (ptrdiff_t) round(mapRow);

                if (intRow < 0 || intRow > (ptrdiff_t) nRow-1) {
                    arrayDst[i] = 0;
                }
                else {
                    arrayDst[i] = arraySrc[intRow];
                }
                break;
            case INTER_LINEAR:
                intRow = (ptrdiff_t) std::floor(mapRow);

                intRow1 = intRow + 1;

                bRow  = (intRow  >= 0) && (intRow  <= (ptrdiff_t) nRow-1);
                bRow1 = (intRow1 >= 0) && (intRow1 <= (ptrdiff_t) nRow-1);

                v0 = 0;     v1 = 0;

                if (bRow )  v0 = arraySrc[intRow ];
                if (bRow1)  v1 = arraySrc[intRow1];

                fracRow = mapRow - (T) intRow;

                arrayDst[i] = (1-fracRow) * v0 + fracRow * v1;
                break;
            case INTER_CUBIC:
                ERROR("transform_translate", "unsupported interpolation mode");
                break;
            default:
                ERROR("transform_translate", "unsupported interpolation mode");
        }
    }
}

// instantiation
template
void transform_translate<float >(const float*  const arraySrc, float*  const arrayDst,
                                 size_t nRow,
                                 float  nRowOff,
                                 eInter inter);
template
void transform_translate<double>(const double* const arraySrc, double* const arrayDst,
                                 size_t nRow,
                                 double nRowOff,
                                 eInter inter);

// 2D
template <typename T>
void transform_translate(const T* const arraySrc, T* const arrayDst,
                         size_t nRow,    size_t nCol,
                         T      nRowOff, T      nColOff,
                         eInter inter)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);

    ptrdiff_t    intRow, intCol, intRow1, intCol1;
    bool         bRow, bCol, bRow1, bCol1;
    T            pixRow, pixCol, mapRow, mapCol;
    T            v00, v01, v10, v11;
    T            fracRow, fracCol;

    //
    // do not use OpenMP here, it slows down the program
    //
    for (size_t i = 0; i < nRow*nCol; i++) {
        // compute original pixel's coordinates
        ind2sub(i, nCol, pixRow, pixCol);

        // compute mapped pixel's coordinates
        mapRow = pixRow - nRowOff;
        mapCol = pixCol - nColOff;

        switch (inter) {
            case INTER_NEAREST:
                intRow = (ptrdiff_t) round(mapRow);
                intCol = (ptrdiff_t) round(mapCol);

                if (intRow < 0 || intRow > (ptrdiff_t) nRow-1 ||
                    intCol < 0 || intCol > (ptrdiff_t) nCol-1) {
                    arrayDst[i] = 0;
                }
                else {
                    arrayDst[i] = arraySrc[intRow*nCol+intCol];
                }
                break;
            case INTER_LINEAR:
                intRow = (ptrdiff_t) std::floor(mapRow);
                intCol = (ptrdiff_t) std::floor(mapCol);

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

                fracRow = mapRow - (T) intRow;
                fracCol = mapCol - (T) intCol;

                arrayDst[i] = (1-fracRow) * (1-fracCol) * v00 +
                              (1-fracRow) * fracCol     * v01 +
                              fracRow     * (1-fracCol) * v10 +
                              fracRow     * fracCol     * v11;
                break;
            case INTER_CUBIC:
                ERROR("transform_translate", "unsupported interpolation mode");
                break;
            default:
                ERROR("transform_translate", "unsupported interpolation mode");
        }
    }
}

// instantiation
template
void transform_translate<float >(const float*  const arraySrc, float*  const arrayDst,
                                 size_t nRow,    size_t nCol,
                                 float  nRowOff, float  nColOff,
                                 eInter inter);
template
void transform_translate<double>(const double* const arraySrc, double* const arrayDst,
                                 size_t nRow,    size_t nCol,
                                 double nRowOff, double nColOff,
                                 eInter inter);

// 3D
template <typename T>
void transform_translate(const T* const arraySrc, T* const arrayDst,
                         size_t nRow,    size_t nCol,    size_t nSec,
                         T      nRowOff, T      nColOff, T      nSecOff,
                         eInter inter)
{
    assert(arraySrc != NULL && arrayDst != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    ptrdiff_t    intRow, intCol, intSec, intRow1, intCol1, intSec1;
    bool         bRow, bCol, bSec, bRow1, bCol1, bSec1;
    T            pixRow, pixCol, pixSec, mapRow, mapCol, mapSec;
    T            v000, v001, v010, v011, v100, v101, v110, v111;
    T            fracRow, fracCol, fracSec;

    //
    // do not use OpenMP here, it slows down the program
    //
    for (size_t i = 0; i < nRow*nCol*nSec; i++) {
        // compute original pixel's coordinates
        ind2sub(i, nCol, nSec, pixRow, pixCol, pixSec);

        // compute mapped pixel's coordinates
        mapRow = pixRow - nRowOff;
        mapCol = pixCol - nColOff;
        mapSec = pixSec - nSecOff;

        switch (inter) {
            case INTER_NEAREST:
                intRow = (ptrdiff_t) round(mapRow);
                intCol = (ptrdiff_t) round(mapCol);
                intSec = (ptrdiff_t) round(mapSec);

                if (intRow < 0 || intRow > (ptrdiff_t) nRow-1 ||
                    intCol < 0 || intCol > (ptrdiff_t) nCol-1 ||
                    intSec < 0 || intSec > (ptrdiff_t) nSec-1) {
                    arrayDst[i] = 0;
                }
                else {
                    arrayDst[i] = arraySrc[intRow*nCol*nSec+intCol*nSec+intSec];
                }
                break;
            case INTER_LINEAR:
                intRow = (ptrdiff_t) std::floor(mapRow);
                intCol = (ptrdiff_t) std::floor(mapCol);
                intSec = (ptrdiff_t) std::floor(mapSec);

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

                fracRow = mapRow - (T) intRow;
                fracCol = mapCol - (T) intCol;
                fracSec = mapSec - (T) intSec;

                arrayDst[i] = (1-fracRow) * (1-fracCol) * (1-fracSec) * v000
                            + (1-fracRow) * (1-fracCol) * fracSec     * v001
                            + (1-fracRow) * fracCol     * (1-fracSec) * v010
                            + (1-fracRow) * fracCol     * fracSec     * v011
                            + fracRow     * (1-fracCol) * (1-fracSec) * v100
                            + fracRow     * (1-fracCol) * fracSec     * v101
                            + fracRow     * fracCol     * (1-fracSec) * v110
                            + fracRow     * fracCol     * fracSec     * v111;
                break;
            case INTER_CUBIC:
                ERROR("transform_translate", "unsupported interpolation mode");
                break;
            default:
                ERROR("transform_translate", "unsupported interpolation mode");
        }
    }
}

// instantiation
template
void transform_translate<float >(const float*  const arraySrc, float*  const arrayDst,
                                 size_t nRow,    size_t nCol,    size_t nSec,
                                 float  nRowOff, float  nColOff, float  nSecOff,
                                 eInter inter);
template
void transform_translate<double>(const double* const arraySrc, double* const arrayDst,
                                 size_t nRow,    size_t nCol,    size_t nSec,
                                 double nRowOff, double nColOff, double nSecOff,
                                 eInter inter);

} // namespace gem
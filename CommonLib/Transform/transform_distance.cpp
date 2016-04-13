/***********************************************************************
 *  File:       transform_distance.cpp
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
 * Distance / Voronoi
 *************************/

// 1D
template <typename T>
void transform_distance(const T*      const arraySrc, size_t nRow,
                              T*      const arrayDst,
                        const size_t* const labelSrc,
                              size_t* const labelDst)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(((labelSrc != NULL) && (labelDst != NULL)) ||
           ((labelSrc == NULL) && (labelDst == NULL)));

    T           distMax = (T) GEM_DISTANCE_MAX;
    size_t      iPar;           // index of rightmost parabola in lower envelope
    size_t*     v = NULL;       // locations of parabolas in lower envelope
    T*          z = NULL;       // locations of boundaries between parabolas
    T           intersectPoint; // intersection point between parabola

    array_new(v, nRow);
    array_new(z, nRow+1);

    iPar = 0;
    v[0] = 0;
    z[0] = -distMax;
    z[1] = +distMax;

    // compute lower envelope
    for (size_t iRow = 1; iRow < nRow; iRow++) {

        intersectPoint  = ((arraySrc[iRow] + pow2((T) iRow)) - (arraySrc[v[iPar]] + pow2((T) v[iPar]))) / (T) (2*iRow-2*v[iPar]);

        while (intersectPoint <= z[iPar]) {
            iPar--;
            intersectPoint  = ((arraySrc[iRow] + pow2((T) iRow)) - (arraySrc[v[iPar]] + pow2((T) v[iPar]))) / (T) (2*iRow-2*v[iPar]);
        }

        iPar++;
        v[iPar] = iRow;
        z[iPar] = intersectPoint;
        z[iPar+1] = +distMax;
    }

    // fill in values of distance transform
    iPar = 0;
    for (size_t iRow = 0; iRow < nRow; iRow++) {

        while (z[iPar+1] < iRow) {
            iPar++;
        }

        arrayDst[iRow] = (T) pow2(iRow-v[iPar]) + arraySrc[v[iPar]];

        if (labelSrc != NULL) {
            labelDst[iRow] = labelSrc[v[iPar]];
        }
    }

    array_delete(v);
    array_delete(z);
}

// instantiation
template
void transform_distance<float >(const float*  const arraySrc, size_t nRow,
                                      float*  const arrayDst,
                                const size_t* const labelSrc,
                                      size_t* const labelDst);
template
void transform_distance<double>(const double* const arraySrc, size_t nRow,
                                      double* const arrayDst,
                                const size_t* const labelSrc,
                                      size_t* const labelDst);

// 2D
template <typename T>
void transform_distance(const T*      const arraySrc, size_t nRow, size_t nCol,
                              T*      const arrayDst,
                        const size_t* const labelSrc,
                              size_t* const labelDst)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
    assert(((labelSrc != NULL) && (labelDst != NULL)) ||
           ((labelSrc == NULL) && (labelDst == NULL)));

    size_t    iRow, iCol;
    T*        arraySrc1D = NULL;
    T*        arrayDst1D = NULL;
    size_t*   labelSrc1D = NULL;
    size_t*   labelDst1D = NULL;

    array_new(arraySrc1D, std::max(nRow,nCol));
    array_new(arrayDst1D, std::max(nRow,nCol));

    if (labelSrc != NULL) {
        array_new(labelSrc1D, std::max(nRow,nCol));
        array_new(labelDst1D, std::max(nRow,nCol));
    }

    // transform along columns
    for (iRow = 0; iRow < nRow; iRow++) {

        for (iCol = 0; iCol < nCol; iCol++) {
            arraySrc1D[iCol] = arraySrc[sub2ind(iRow,iCol,nCol)];
        }

        if (labelSrc != NULL) {
            for (iCol = 0; iCol < nCol; iCol++) {
                labelSrc1D[iCol] = labelSrc[sub2ind(iRow,iCol,nCol)];
            }
        }

        transform_distance(arraySrc1D, nCol,
                           arrayDst1D,
                           labelSrc1D,
                           labelDst1D);

        for (iCol = 0; iCol < nCol; iCol++) {
            arrayDst[sub2ind(iRow,iCol,nCol)] = arrayDst1D[iCol];
        }

        if (labelSrc != NULL) {
            for (iCol = 0; iCol < nCol; iCol++) {
                labelDst[sub2ind(iRow,iCol,nCol)] = labelDst1D[iCol];
            }
        }
    }

    // transform along rows
    for (iCol = 0; iCol < nCol; iCol++) {

        for (iRow = 0; iRow < nRow; iRow++) {
            arraySrc1D[iRow] = arrayDst[sub2ind(iRow,iCol,nCol)];
        }

        if (labelSrc != NULL) {
            for (iRow = 0; iRow < nRow; iRow++) {
                labelSrc1D[iRow] = labelDst[sub2ind(iRow,iCol,nCol)];
            }
        }

        transform_distance(arraySrc1D, nRow,
                           arrayDst1D,
                           labelSrc1D,
                           labelDst1D);

        for (iRow = 0; iRow < nRow; iRow++) {
            arrayDst[sub2ind(iRow,iCol,nCol)] = arrayDst1D[iRow];
        }

        if (labelSrc != NULL) {
            for (iRow = 0; iRow < nRow; iRow++) {
                labelDst[sub2ind(iRow,iCol,nCol)] = labelDst1D[iRow];
            }
        }
    }

    array_delete(arraySrc1D);
    array_delete(arrayDst1D);

    if (labelSrc != NULL) {
        array_delete(labelSrc1D);
        array_delete(labelDst1D);
    }
}

// instantiation
template
void transform_distance<float >(const float*  const arraySrc, size_t nRow, size_t nCol,
                                      float*  const arrayDst,
                                const size_t* const labelSrc,
                                      size_t* const labelDst);
template
void transform_distance<double>(const double* const arraySrc, size_t nRow, size_t nCol,
                                      double* const arrayDst,
                                const size_t* const labelSrc,
                                      size_t* const labelDst);

// 3D
template <typename T>
void transform_distance(const T*      const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                              T*      const arrayDst,
                        const size_t* const labelSrc,
                              size_t* const labelDst)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
    assert(nSec > 0);
    assert(((labelSrc != NULL) && (labelDst != NULL)) ||
           ((labelSrc == NULL) && (labelDst == NULL)));

    size_t    iRow, iCol, iSec;
    T*        arraySrc1D = NULL;
    T*        arrayDst1D = NULL;
    size_t*   labelSrc1D = NULL;
    size_t*   labelDst1D = NULL;

    array_new(arraySrc1D, std::max(std::max(nRow,nCol),nSec));
    array_new(arrayDst1D, std::max(std::max(nRow,nCol),nSec));

    if (labelSrc != NULL) {
        array_new(labelSrc1D, std::max(std::max(nRow,nCol),nSec));
        array_new(labelDst1D, std::max(std::max(nRow,nCol),nSec));
    }

    // transform along sections
    for (iRow = 0; iRow < nRow; iRow++) {
        for (iCol = 0; iCol < nCol; iCol++) {

            for (iSec = 0; iSec < nSec; iSec++) {
                arraySrc1D[iSec] = arraySrc[sub2ind(iRow,iCol,iSec,nCol,nSec)];
            }

            if (labelSrc != NULL) {
                for (iSec = 0; iSec < nSec; iSec++) {
                    labelSrc1D[iSec] = labelSrc[sub2ind(iRow,iCol,iSec,nCol,nSec)];
                }
            }

            transform_distance(arraySrc1D, nSec,
                               arrayDst1D,
                               labelSrc1D,
                               labelDst1D);

            for (iSec = 0; iSec < nSec; iSec++) {
                arrayDst[sub2ind(iRow,iCol,iSec,nCol,nSec)] = arrayDst1D[iSec];
            }

            if (labelSrc != NULL) {
                for (iSec = 0; iSec < nSec; iSec++) {
                    labelDst[sub2ind(iRow,iCol,iSec,nCol,nSec)] = labelDst1D[iSec];
                }
            }
        }
    }

    // transform along columns
    for (iRow = 0; iRow < nRow; iRow++) {
        for (iSec = 0; iSec < nSec; iSec++) {

            for (iCol = 0; iCol < nCol; iCol++) {
                arraySrc1D[iCol] = arrayDst[sub2ind(iRow,iCol,iSec,nCol,nSec)];
            }

            if (labelSrc != NULL) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    labelSrc1D[iCol] = labelDst[sub2ind(iRow,iCol,iSec,nCol,nSec)];
                }
            }

            transform_distance(arraySrc1D, nCol,
                               arrayDst1D,
                               labelSrc1D,
                               labelDst1D);

            for (iCol = 0; iCol < nCol; iCol++) {
                arrayDst[sub2ind(iRow,iCol,iSec,nCol,nSec)] = arrayDst1D[iCol];
            }

            if (labelSrc != NULL) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    labelDst[sub2ind(iRow,iCol,iSec,nCol,nSec)] = labelDst1D[iCol];
                }
            }
        }
    }

    // transform along rows
    for (iCol = 0; iCol < nCol; iCol++) {
        for (iSec = 0; iSec < nSec; iSec++) {

            for (iRow = 0; iRow < nRow; iRow++) {
                arraySrc1D[iRow] = arrayDst[sub2ind(iRow,iCol,iSec,nCol,nSec)];
            }

            if (labelSrc != NULL) {
                for (iRow = 0; iRow < nRow; iRow++) {
                    labelSrc1D[iRow] = labelDst[sub2ind(iRow,iCol,iSec,nCol,nSec)];
                }
            }

            transform_distance(arraySrc1D, nCol,
                               arrayDst1D,
                               labelSrc1D,
                               labelDst1D);

            for (iRow = 0; iRow < nRow; iRow++) {
                arrayDst[sub2ind(iRow,iCol,iSec,nCol,nSec)] = arrayDst1D[iRow];
            }

            if (labelSrc != NULL) {
                for (iRow = 0; iRow < nRow; iRow++) {
                    labelDst[sub2ind(iRow,iCol,iSec,nCol,nSec)] = labelDst1D[iRow];
                }
            }
        }
    }

    array_delete(arraySrc1D);
    array_delete(arrayDst1D);

    if (labelSrc != NULL) {
        array_delete(labelSrc1D);
        array_delete(labelDst1D);
    }
}

// instantiation
template
void transform_distance<float >(const float*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      float*  const arrayDst,
                                const size_t* const labelSrc,
                                      size_t* const labelDst);
template
void transform_distance<double>(const double* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                      double* const arrayDst,
                                const size_t* const labelSrc,
                                      size_t* const labelDst);

} // namespace gem

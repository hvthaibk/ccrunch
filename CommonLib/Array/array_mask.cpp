/***********************************************************************
 *  File:       array_mask.cpp
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
 * Mask
 ****************************************/

// 1D copy
#ifndef NDEBUG
template <typename T>
void array_mask_copy(const T* const array,  size_t    nRow,
                     const T* const mask,   size_t    nRowMsk,
                           T* const output, ptrdiff_t nRowOff)
#else
template <typename T>
void array_mask_copy(const T* const array,  size_t    ,
                     const T* const mask,   size_t    nRowMsk,
                           T* const output, ptrdiff_t nRowOff)
#endif
{
    assert(array  != NULL);
    assert(mask   != NULL);
    assert(output != NULL);
    assert(nRow > 0 && nRowMsk > 0);

    ptrdiff_t    nRowMskSigned = (ptrdiff_t) nRowMsk;
    ptrdiff_t    iRowMsk;

    #pragma omp parallel for
    for (iRowMsk = 0; iRowMsk < nRowMskSigned; iRowMsk++) {
        if (mask[iRowMsk] > 0) {
            assert(iRowMsk + nRowOff >= 0);
            assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);

            output[iRowMsk] = array[iRowMsk+nRowOff];
        }
        else {
            output[iRowMsk] = 0;
        }
    }
}

// instantiation
template
void array_mask_copy<float >(const float*  const array,  size_t    nRow,
                             const float*  const mask,   size_t    nRowMsk,
                                   float*  const output, ptrdiff_t nRowOff);
template
void array_mask_copy<double>(const double* const array,  size_t    nRow,
                             const double* const mask,   size_t    nRowMsk,
                                   double* const output, ptrdiff_t nRowOff);

// 2D copy
#ifndef NDEBUG
template <typename T>
void array_mask_copy(const T* const array,  size_t    nRow,    size_t    nCol,
                     const T* const mask,   size_t    nRowMsk, size_t    nColMsk,
                           T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff)
#else
template <typename T>
void array_mask_copy(const T* const array,  size_t        ,    size_t    nCol,
                     const T* const mask,   size_t    nRowMsk, size_t    nColMsk,
                           T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff)
#endif
{
    assert(array  != NULL);
    assert(mask   != NULL);
    assert(output != NULL);
    assert(nRow > 0 && nRowMsk > 0);
    assert(nCol > 0 && nColMsk > 0);

    ptrdiff_t    nRowMskSigned = (ptrdiff_t) nRowMsk;
    ptrdiff_t    nColMskSigned = (ptrdiff_t) nColMsk;
    ptrdiff_t    iRowMsk, iColMsk;
    ptrdiff_t    iSrc1;
    ptrdiff_t    iMsk1;
    ptrdiff_t    nOff = nRowOff*nCol + nColOff;

    #pragma omp parallel for private(iRowMsk,iSrc1,iMsk1,iColMsk)
    for (iRowMsk = 0; iRowMsk < nRowMskSigned; iRowMsk++) {
        iSrc1 = iRowMsk*nCol + nOff;
        iMsk1 = iRowMsk*nColMsk;

        for (iColMsk = 0; iColMsk < nColMskSigned; iColMsk++) {
            if (mask[iMsk1+iColMsk] > 0) {
                assert(iRowMsk + nRowOff >= 0);
                assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);
                assert(iColMsk + nColOff >= 0);
                assert(iColMsk + nColOff <  (ptrdiff_t) nCol);

                output[iMsk1+iColMsk] = array[iSrc1+iColMsk];
            }
            else {
                output[iMsk1+iColMsk] = 0;
            }
        }
    }
}

// instantiation
template
void array_mask_copy<float >(const float*  const array,  size_t    nRow,    size_t    nCol,
                             const float*  const mask,   size_t    nRowMsk, size_t    nColMsk,
                                   float*  const output, ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void array_mask_copy<double>(const double* const array,  size_t    nRow,    size_t    nCol,
                             const double* const mask,   size_t    nRowMsk, size_t    nColMsk,
                                   double* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D copy
#ifndef NDEBUG
template <typename T>
void array_mask_copy(const T* const array,  size_t    nRow,    size_t    nCol,    size_t    nSec,
                     const T* const mask,   size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                           T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff)
#else
template <typename T>
void array_mask_copy(const T* const array,  size_t        ,    size_t    nCol,    size_t    nSec,
                     const T* const mask,   size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                           T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff)
#endif
{
    assert(array  != NULL);
    assert(mask   != NULL);
    assert(output != NULL);
    assert(nRow > 0 && nRowMsk > 0);
    assert(nCol > 0 && nColMsk > 0);
    assert(nSec > 0 && nSecMsk > 0);

    ptrdiff_t    nRowMskSigned = (ptrdiff_t) nRowMsk;
    ptrdiff_t    nColMskSigned = (ptrdiff_t) nColMsk;
    ptrdiff_t    nSecMskSigned = (ptrdiff_t) nSecMsk;
    ptrdiff_t    iRowMsk, iColMsk, iSecMsk;
    ptrdiff_t    iSrc1, iSrc2;
    ptrdiff_t    iMsk1, iMsk2;
    ptrdiff_t    nColSec = nCol*nSec;
    ptrdiff_t    nColSecMsk = nColMsk*nSecMsk;
    ptrdiff_t    nOff = (nRowOff*nCol + nColOff)*nSec + nSecOff;

    #pragma omp parallel for private(iRowMsk,iSrc1,iMsk1,iColMsk,iSrc2,iMsk2,iSecMsk)
    for (iRowMsk = 0; iRowMsk < nRowMskSigned; iRowMsk++) {
        iSrc1 = iRowMsk*nColSec + nOff;
        iMsk1 = iRowMsk*nColSecMsk;

        for (iColMsk = 0; iColMsk < nColMskSigned; iColMsk++) {
            iSrc2 = iSrc1 + iColMsk*nSec;
            iMsk2 = iMsk1 + iColMsk*nSecMsk;

            for (iSecMsk = 0; iSecMsk < nSecMskSigned; iSecMsk++) {
                if (mask[iMsk2+iSecMsk] > 0) {
                    assert(iRowMsk + nRowOff >= 0);
                    assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);
                    assert(iColMsk + nColOff >= 0);
                    assert(iColMsk + nColOff <  (ptrdiff_t) nCol);
                    assert(iSecMsk + nSecOff >= 0);
                    assert(iSecMsk + nSecOff <  (ptrdiff_t) nSec);

                    output[iMsk2+iSecMsk] = array[iSrc2+iSecMsk];
                }
                else {
                    output[iMsk2+iSecMsk] = 0;
                }
            }
        }
    }
}

// instantiation
template
void array_mask_copy<float >(const float*  const array,  size_t    nRow,    size_t    nCol,    size_t    nSec,
                             const float*  const mask,   size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                   float*  const output, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void array_mask_copy<double>(const double* const array,  size_t    nRow,    size_t    nCol,    size_t    nSec,
                             const double* const mask,   size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                   double* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

// 1D replace
//#ifndef NDEBUG
template <typename T>
void array_mask_replace(      T* const array, size_t    nRow,
                        const T* const mask,  size_t    nRowMsk,
                              T        value, ptrdiff_t nRowOff,
                              bool     maskInside)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0 && nRowMsk > 0);

    ptrdiff_t    nRowMskSigned = (ptrdiff_t) nRowMsk;
    ptrdiff_t    iRowMsk;

    #pragma omp parallel for
    for (iRowMsk = 0; iRowMsk < nRowMskSigned; iRowMsk++) {
        if (mask[iRowMsk] > 0) {
            if (maskInside) {
                assert(iRowMsk + nRowOff >= 0);
                assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);

                array[iRowMsk+nRowOff] = value;
            }
            else {
                if (iRowMsk + nRowOff >= 0 &&
                    iRowMsk + nRowOff < (ptrdiff_t) nRow) {
                    array[iRowMsk+nRowOff] = value;
                }
            }
        }
    }
}

// instantiation
template
void array_mask_replace<float >(      float*  const array, size_t    nRow,
                                const float*  const mask,  size_t    nRowMsk,
                                      float         value, ptrdiff_t nRowOff,
                                      bool          maskInside);
template
void array_mask_replace<double>(      double* const array, size_t    nRow,
                                const double* const mask,  size_t    nRowMsk,
                                      double        value, ptrdiff_t nRowOff,
                                      bool          maskInside);

// 2D replace
//#ifndef NDEBUG
template <typename T>
void array_mask_replace(      T* const array, size_t    nRow,    size_t    nCol,
                        const T* const mask,  size_t    nRowMsk, size_t    nColMsk,
                              T        value, ptrdiff_t nRowOff, ptrdiff_t nColOff,
                              bool     maskInside)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0 && nRowMsk > 0);
    assert(nCol > 0 && nColMsk > 0);

    ptrdiff_t    nRowMskSigned = (ptrdiff_t) nRowMsk;
    ptrdiff_t    nColMskSigned = (ptrdiff_t) nColMsk;
    ptrdiff_t    iRowMsk, iColMsk;
    ptrdiff_t    iSrc1;
    ptrdiff_t    iMsk1;
    ptrdiff_t    nOff = nRowOff*nCol + nColOff;

    #pragma omp parallel for private(iRowMsk,iSrc1,iMsk1,iColMsk)
    for (iRowMsk = 0; iRowMsk < nRowMskSigned; iRowMsk++) {
        iSrc1 = iRowMsk*nCol + nOff;
        iMsk1 = iRowMsk*nColMsk;

        for (iColMsk = 0; iColMsk < nColMskSigned; iColMsk++) {
            if (mask[iMsk1+iColMsk] > 0) {
                if (maskInside) {
                    assert(iRowMsk + nRowOff >= 0);
                    assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);
                    assert(iColMsk + nColOff >= 0);
                    assert(iColMsk + nColOff <  (ptrdiff_t) nCol);

                    array[iSrc1+iColMsk] = value;
                }
                else {
                    if (iRowMsk + nRowOff >= 0 &&
                        iRowMsk + nRowOff < (ptrdiff_t) nRow &&
                        iColMsk + nColOff >= 0 &&
                        iColMsk + nColOff < (ptrdiff_t) nCol) {
                        array[iSrc1+iColMsk] = value;
                    }
                }
            }
        }
    }
}

// instantiation
template
void array_mask_replace<float >(      float*  const array, size_t    nRow,    size_t    nCol,
                                const float*  const mask,  size_t    nRowMsk, size_t    nColMsk,
                                      float         value, ptrdiff_t nRowOff, ptrdiff_t nColOff,
                                      bool          maskInside);
template
void array_mask_replace<double>(      double* const array, size_t    nRow,    size_t    nCol,
                                const double* const mask,  size_t    nRowMsk, size_t    nColMsk,
                                      double        value, ptrdiff_t nRowOff, ptrdiff_t nColOff,
                                      bool          maskInside);

// 3D replace
//#ifndef NDEBUG
template <typename T>
void array_mask_replace(      T* const array, size_t    nRow,    size_t    nCol,    size_t    nSec,
                        const T* const mask,  size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                              T        value, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff,
                              bool     maskInside)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0 && nRowMsk > 0);
    assert(nCol > 0 && nColMsk > 0);
    assert(nSec > 0 && nSecMsk > 0);

    ptrdiff_t    nRowMskSigned = (ptrdiff_t) nRowMsk;
    ptrdiff_t    nColMskSigned = (ptrdiff_t) nColMsk;
    ptrdiff_t    nSecMskSigned = (ptrdiff_t) nSecMsk;
    ptrdiff_t    iRowMsk, iColMsk, iSecMsk;
    ptrdiff_t    iSrc1, iSrc2;
    ptrdiff_t    iMsk1, iMsk2;
    ptrdiff_t    nColSec = nCol*nSec;
    ptrdiff_t    nColSecMsk = nColMsk*nSecMsk;
    ptrdiff_t    nOff = (nRowOff*nCol + nColOff)*nSec + nSecOff;

    #pragma omp parallel for private(iRowMsk,iSrc1,iMsk1,iColMsk,iSrc2,iMsk2,iSecMsk)
    for (iRowMsk = 0; iRowMsk < nRowMskSigned; iRowMsk++) {
        iSrc1 = iRowMsk*nColSec + nOff;
        iMsk1 = iRowMsk*nColSecMsk;

        for (iColMsk = 0; iColMsk < nColMskSigned; iColMsk++) {
            iSrc2 = iSrc1 + iColMsk*nSec;
            iMsk2 = iMsk1 + iColMsk*nSecMsk;

            for (iSecMsk = 0; iSecMsk < nSecMskSigned; iSecMsk++) {
                if (mask[iMsk2+iSecMsk] > 0) {
                    if (maskInside) {
                        assert(iRowMsk + nRowOff >= 0);
                        assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);
                        assert(iColMsk + nColOff >= 0);
                        assert(iColMsk + nColOff <  (ptrdiff_t) nCol);
                        assert(iSecMsk + nSecOff >= 0);
                        assert(iSecMsk + nSecOff <  (ptrdiff_t) nSec);

                        array[iSrc2+iSecMsk] = value;
                    }
                    else {
                        if (iRowMsk + nRowOff >= 0 &&
                            iRowMsk + nRowOff < (ptrdiff_t) nRow &&
                            iColMsk + nColOff >= 0 &&
                            iColMsk + nColOff < (ptrdiff_t) nCol &&
                            iSecMsk + nSecOff >= 0 &&
                            iSecMsk + nSecOff < (ptrdiff_t) nSec) {
                            array[iSrc2+iSecMsk] = value;
                        }
                    }
                }
            }
        }
    }
}

// instantiation
template
void array_mask_replace<float >(      float*  const array, size_t    nRow,    size_t    nCol,    size_t    nSec,
                                const float*  const mask,  size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                      float         value, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff,
                                      bool          maskInside);
template
void array_mask_replace<double>(      double* const array, size_t    nRow,    size_t    nCol,    size_t    nSec,
                                const double* const mask,  size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                      double        value, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff,
                                      bool          maskInside);

// sum
template <typename T>
T array_mask_sum(const T* const array, const T* const mask, size_t nRow)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0);

    T        sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (mask[iRow] > 0) {
            sum += array[iRow];
        }
    }

    return sum;
}

// instantiation
template
float  array_mask_sum<float >(const float*  const array, const float*  const mask, size_t nRow);
template
double array_mask_sum<double>(const double* const array, const double* const mask, size_t nRow);

// substract
template <typename T>
void array_mask_sub(T* const arrayDst, const T* const arraySrc, const T* const mask, T value, size_t nRow)
{
    assert(arrayDst != NULL);
    assert(arraySrc != NULL);
    assert(mask     != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (mask[iRow] > 0) {
            arrayDst[iRow] = arraySrc[iRow] - value;
        }
    }
}

// instantiation
template
void array_mask_sub(float*  const arrayDst, const float*  const arraySrc, const float*  const mask, float  value, size_t nRow);
template
void array_mask_sub(double* const arrayDst, const double* const arraySrc, const double* const mask, double value, size_t nRow);

// substract - square - sum
template <typename T>
T array_mask_subsqrsum(const T* const array, T value, const T* const mask, size_t nRow)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0);

    T        sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (mask[iRow] > 0) {
            sum += pow2(array[iRow]-value);
        }
    }

    return sum;
}

// instantiation
template
float  array_mask_subsqrsum(const float*  const array, float  value, const float*  const mask, size_t nRow);
template
double array_mask_subsqrsum(const double* const array, double value, const double* const mask, size_t nRow);

// substract - divide - multiply
template <typename T>
void array_mask_subdivmul(T* const array, T valsub, T valdiv, const T* const mask, size_t nRow)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] = (array[iRow] - valsub) / valdiv * mask[iRow];
    }
}

// instantiation
template
void array_mask_subdivmul(float*  const array, float  valsub, float  valdiv, const float*  const mask, size_t nRow);
template
void array_mask_subdivmul(double* const array, double valsub, double valdiv, const double* const mask, size_t nRow);

} // namespace gem

/***********************************************************************
 *  File:       filter_dog.cpp
 *
 *  Purpose:    Implementation of filtering functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "filter.hpp"

namespace gem {

/**************************
 * Difference of Gaussian
 *************************/

// 1D
template <typename T>
void filter_dog(T* const array, size_t nRow, T sigma1, T sigma2)
{
    assert(array != NULL);
    assert(nRow > 0);
    assert(sigma1 > 0 && sigma2 > 0);

    T        fPre1 = 1 / (sigma1 * (T) sqrt(2*M_PI));
    T        fExp1 = -1 / (2 * pow2(sigma1));
    T        fPre2 = 1 / (sigma2 * (T) sqrt(2*M_PI));
    T        fExp2 = -1 / (2 * pow2(sigma2));
    T        cRow = ((T) nRow - 1) / 2;
    T        *array1 = NULL, *array2 = NULL;
    T        distToCenterSqr;
    size_t   iRow;

    array_new(array1, nRow);
    array_new(array2, nRow);

    #pragma omp parallel for private(iRow,distToCenterSqr)
    for (iRow = 0; iRow < nRow; iRow++) {
        distToCenterSqr = distSqr((T)iRow,cRow);

        array1[iRow] = fPre1 * std::exp(fExp1 * distToCenterSqr);
        array2[iRow] = fPre2 * std::exp(fExp2 * distToCenterSqr);
    }

    T    arraySum;
    arraySum = array_reduce_sum(array1, nRow);
    array_math_div(array1, arraySum, nRow);
    arraySum = array_reduce_sum(array2, nRow);
    array_math_div(array2, arraySum, nRow);

    array_math_sub(array, array1, array2, nRow);

    array_delete(array1);
    array_delete(array2);
}

// instantiation
template
void filter_dog<float >(float*  const array, size_t nRow, float  sigma1, float  sigma2);
template
void filter_dog<double>(double* const array, size_t nRow, double sigma1, double sigma2);

// 2D
template <typename T>
void filter_dog(T* const array, size_t nRow, size_t nCol, T sigma1, T sigma2)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0);
    assert(sigma1 > 0 && sigma2 > 0);

    T        fPre1 = 1 / (sigma1 * (T) sqrt(2*M_PI));
    T        fExp1 = -1 / (2 * pow2(sigma1));
    T        fPre2 = 1 / (sigma2 * (T) sqrt(2*M_PI));
    T        fExp2 = -1 / (2 * pow2(sigma2));
    T        cRow = ((T) nRow - 1) / 2;
    T        cCol = ((T) nCol - 1) / 2;
    T        *array1 = NULL, *array2 = NULL;
    T        dist1, distToCenterSqr, distSqrMax;
    size_t   indx1;
    size_t   iRow, iCol;

    distSqrMax = std::max(pow2(cRow),pow2(cCol));

    array_new(array1, nRow*nCol);
    array_new(array2, nRow*nCol);

    #pragma omp parallel for private(iRow,dist1,indx1,iCol,distToCenterSqr)
    for (iRow = 0; iRow < nRow; iRow++) {
        dist1 = distSqr((T)iRow,cRow);
        indx1 = iRow*nCol;

        for (iCol = 0; iCol < nCol; iCol++) {
            distToCenterSqr = dist1 + distSqr((T)iCol,cCol);

            if (distToCenterSqr <= distSqrMax) {
                array1[indx1+iCol] = fPre1 * std::exp(fExp1 * distToCenterSqr);
                array2[indx1+iCol] = fPre2 * std::exp(fExp2 * distToCenterSqr);
            }
            else {
                array1[indx1+iCol] = 0;
                array2[indx1+iCol] = 0;
            }
        }
    }

    T    arraySum;
    arraySum = array_reduce_sum(array1, nRow*nCol);
    array_math_div(array1, arraySum, nRow*nCol);
    arraySum = array_reduce_sum(array2, nRow*nCol);
    array_math_div(array2, arraySum, nRow*nCol);

    array_math_sub(array, array1, array2, nRow*nCol);

    array_delete(array1);
    array_delete(array2);
}

// instantiation
template
void filter_dog<float >(float*  const array, size_t nRow, size_t nCol, float  sigma1, float  sigma2);
template
void filter_dog<double>(double* const array, size_t nRow, size_t nCol, double sigma1, double sigma2);

// 3D
template <typename T>
void filter_dog(T* const array, size_t nRow, size_t nCol, size_t nSec, T sigma1, T sigma2)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);
    assert(sigma1 > 0 && sigma2 > 0);

    T        fPre1 = 1 / (sigma1 * (T) sqrt(2*M_PI));
    T        fExp1 = -1 / (2 * pow2(sigma1));
    T        fPre2 = 1 / (sigma2 * (T) sqrt(2*M_PI));
    T        fExp2 = -1 / (2 * pow2(sigma2));
    T        cRow = ((T) nRow - 1) / 2;
    T        cCol = ((T) nCol - 1) / 2;
    T        cSec = ((T) nSec - 1) / 2;
    T        *array1 = NULL, *array2 = NULL;
    T        dist1, dist2, distToCenterSqr, distSqrMax;
    size_t   indx1, indx2;
    size_t   iRow, iCol, iSec;

    distSqrMax = std::max(std::max(pow2(cRow),pow2(cCol)),pow2(cSec));

    array_new(array1, nRow*nCol*nSec);
    array_new(array2, nRow*nCol*nSec);

    #pragma omp parallel for private(iRow,dist1,indx1,iCol,dist2,indx2,iSec,distToCenterSqr)
    for (iRow = 0; iRow < nRow; iRow++) {
        dist1 = distSqr((T)iRow,cRow);
        indx1 = iRow*nCol*nSec;

        for (iCol = 0; iCol < nCol; iCol++) {
            dist2 = dist1 + distSqr((T)iCol,cCol);
            indx2 = indx1 + iCol*nSec;

            for (iSec = 0; iSec < nSec; iSec++) {
                distToCenterSqr = dist2 + distSqr((T)iSec,cSec);;

                if (distToCenterSqr <= distSqrMax) {
                    array1[indx2+iSec] = fPre1 * std::exp(fExp1 * distToCenterSqr);
                    array2[indx2+iSec] = fPre2 * std::exp(fExp2 * distToCenterSqr);
                }
                else {
                    array1[indx2+iSec] = 0;
                    array2[indx2+iSec] = 0;
                }
            }
        }
    }

    T    arraySum;
    arraySum = array_reduce_sum(array1, nRow*nCol*nSec);
    array_math_div(array1, arraySum, nRow*nCol*nSec);
    arraySum = array_reduce_sum(array2, nRow*nCol*nSec);
    array_math_div(array2, arraySum, nRow*nCol*nSec);

    array_math_sub(array, array1, array2, nRow*nCol*nSec);

    array_delete(array1);
    array_delete(array2);
}

// instantiation
template
void filter_dog<float >(float*  const array, size_t nRow, size_t nCol, size_t nSec, float  sigma1, float  sigma2);
template
void filter_dog<double>(double* const array, size_t nRow, size_t nCol, size_t nSec, double sigma1, double sigma2);

} // namespace gem

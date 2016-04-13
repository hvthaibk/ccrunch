/***********************************************************************
 *  File:       filter_gaussian.cpp
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
 * Gaussian
 *************************/

// 1D
template <typename T>
void filter_gaussian(T* const array, size_t nRow, T sigma)
{
    assert(array != NULL);
    assert(nRow > 0);
    assert(sigma > 0);

    T        fExp = -1 / (2 * pow2(sigma));
    T        cRow = ((T) nRow - 1) / 2;
    T        distToCenterSqr;
    size_t   iRow;

    #pragma omp parallel for private(iRow,distToCenterSqr)
    for (iRow = 0; iRow < nRow; iRow++) {
        distToCenterSqr = distSqr((T)iRow,cRow);

        array[iRow] = std::exp(fExp * distToCenterSqr);
    }

    T    arraySum = array_reduce_sum(array, nRow);
    array_math_div(array, arraySum, nRow);
}

// instantiation
template
void filter_gaussian<float >(float*  const array, size_t nRow, float  sigma);
template
void filter_gaussian<double>(double* const array, size_t nRow, double sigma);

// 2D
template <typename T>
void filter_gaussian(T* const array, size_t nRow, size_t nCol, T sigma)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0);
    assert(sigma > 0);

    T        fExp = -1 / (2 * pow2(sigma));
    T        cRow = ((T) nRow - 1) / 2;
    T        cCol = ((T) nCol - 1) / 2;
    T        dist1, distToCenterSqr, distSqrMax;
    size_t   indx1;
    size_t   iRow, iCol;

    distSqrMax = std::max(pow2(cRow),pow2(cCol));

    #pragma omp parallel for private(iRow,dist1,indx1,iCol,distToCenterSqr)
    for (iRow = 0; iRow < nRow; iRow++) {
        dist1 = distSqr((T)iRow,cRow);
        indx1 = iRow*nCol;

        for (iCol = 0; iCol < nCol; iCol++) {
            distToCenterSqr = dist1 + distSqr((T)iCol,cCol);

            if (distToCenterSqr <= distSqrMax) {
                array[indx1+iCol] = std::exp(fExp * distToCenterSqr);
            }
            else {
                array[indx1+iCol] = 0;
            }
        }
    }

    T    arraySum = array_reduce_sum(array, nRow*nCol);
    array_math_div(array, arraySum, nRow*nCol);
}

// instantiation
template
void filter_gaussian<float >(float*  const array, size_t nRow, size_t nCol, float  sigma);
template
void filter_gaussian<double>(double* const array, size_t nRow, size_t nCol, double sigma);

// 3D
template <typename T>
void filter_gaussian(T* const array, size_t nRow, size_t nCol, size_t nSec, T sigma)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);
    assert(sigma > 0);

    T        fExp = -1 / (2 * pow2(sigma));
    T        cRow = ((T) nRow - 1) / 2;
    T        cCol = ((T) nCol - 1) / 2;
    T        cSec = ((T) nSec - 1) / 2;
    T        dist1, dist2, distToCenterSqr, distSqrMax;
    size_t   indx1, indx2;
    size_t   iRow, iCol, iSec;

    distSqrMax = std::max(std::max(pow2(cRow),pow2(cCol)),pow2(cSec));

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
                    array[indx2+iSec] = std::exp(fExp * distToCenterSqr);
                }
                else {
                    array[indx2+iSec] = 0;
                }
            }
        }
    }

    T    arraySum = array_reduce_sum(array, nRow*nCol*nSec);
    array_math_div(array, arraySum, nRow*nCol*nSec);
}

// instantiation
template
void filter_gaussian<float >(float*  const array, size_t nRow, size_t nCol, size_t nSec, float  sigma);
template
void filter_gaussian<double>(double* const array, size_t nRow, size_t nCol, size_t nSec, double sigma);

} // namespace gem

/***********************************************************************
 *  File:       filter_log.cpp
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
 * Laplacian of Gaussian
 *************************/

// 1D
template <typename T>
void filter_log(T* const array, size_t nRow, T sigma)
{
    assert(array != NULL);
    assert(nRow > 0);
    assert(sigma > 0);

    T        fExp = -1 / (2 * pow2(sigma));
    T        fAdd = -pow2(sigma);
    T        cRow = ((T) nRow - 1) / 2;
    T        expTmp, expSum = 0;
    T        distToCenterSqr;
    size_t   iRow;

    #pragma omp parallel for reduction(+:expSum) private(iRow,expTmp,distToCenterSqr)
    for (iRow = 0; iRow < nRow; iRow++) {
        distToCenterSqr = distSqr((T)iRow,cRow);
        expTmp = std::exp(fExp * distToCenterSqr);
        expSum += expTmp;

        array[iRow] = expTmp * (fAdd + distToCenterSqr);
    }

    T    arrayMean;
    array_math_mul(array, 1/(pow4(sigma)*expSum), nRow);
    arrayMean = array_reduce_mean(array, nRow);
    array_math_sub(array, arrayMean, nRow);
}

// instantiation
template
void filter_log<float >(float*  const array, size_t nRow, float  sigma);
template
void filter_log<double>(double* const array, size_t nRow, double sigma);

// 2D
template <typename T>
void filter_log(T* const array, size_t nRow, size_t nCol, T sigma)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0);
    assert(sigma > 0);

    T        fExp = -1 / (2 * pow2(sigma));
    T        fAdd = -2 * pow2(sigma);
    T        cRow = ((T) nRow - 1) / 2;
    T        cCol = ((T) nCol - 1) / 2;
    T        expTmp, expSum = 0;
    T        dist1, distToCenterSqr;
    size_t   indx1;
    size_t   iRow, iCol;

    #pragma omp parallel for reduction(+:expSum) private(iRow,dist1,indx1,iCol,distToCenterSqr,expTmp)
    for (iRow = 0; iRow < nRow; iRow++) {
        dist1 = distSqr((T)iRow,cRow);
        indx1 = iRow*nCol;

        for (iCol = 0; iCol < nCol; iCol++) {
            distToCenterSqr = dist1 + distSqr((T)iCol,cCol);
            expTmp = std::exp(fExp * distToCenterSqr);
            expSum += expTmp;

            array[indx1+iCol] = expTmp * (fAdd + distToCenterSqr);
        }
    }

    T    arrayMean;
    array_math_mul(array, 1/(pow4(sigma)*expSum), nRow*nCol);
    arrayMean = array_reduce_mean(array, nRow*nCol);
    array_math_sub(array, arrayMean, nRow*nCol);
}

// instantiation
template
void filter_log<float >(float*  const array, size_t nRow, size_t nCol, float  sigma);
template
void filter_log<double>(double* const array, size_t nRow, size_t nCol, double sigma);

// 3D
template <typename T>
void filter_log(T* const array, size_t nRow, size_t nCol, size_t nSec, T sigma)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);
    assert(sigma > 0);

    T        fExp = -1 / (2 * pow2(sigma));
    T        fAdd = -3 * pow2(sigma);
    T        cRow = ((T) nRow - 1) / 2;
    T        cCol = ((T) nCol - 1) / 2;
    T        cSec = ((T) nSec - 1) / 2;
    T        expTmp, expSum = 0;
    T        dist1, dist2, distToCenterSqr;
    size_t   indx1, indx2;
    size_t   iRow, iCol, iSec;

    #pragma omp parallel for reduction(+:expSum) private(iRow,dist1,indx1,iCol,dist2,indx2,iSec,distToCenterSqr,expTmp)
    for (iRow = 0; iRow < nRow; iRow++) {
        dist1 = distSqr((T)iRow,cRow);
        indx1 = iRow*nCol*nSec;

        for (iCol = 0; iCol < nCol; iCol++) {
            dist2 = dist1 + distSqr((T)iCol,cCol);
            indx2 = indx1 + iCol*nSec;

            for (iSec = 0; iSec < nSec; iSec++) {
                distToCenterSqr = dist2 + distSqr((T)iSec,cSec);
                expTmp = std::exp(fExp * distToCenterSqr);
                expSum += expTmp;

                array[indx2+iSec] = expTmp * (fAdd + distToCenterSqr);
            }
        }
    }

    T    arrayMean;
    array_math_mul(array, 1/(pow4(sigma)*expSum), nRow*nCol*nSec);
    arrayMean = array_reduce_mean(array, nRow*nCol*nSec);
    array_math_sub(array, arrayMean, nRow*nCol*nSec);
}

// instantiation
template
void filter_log<float >(float*  const array, size_t nRow, size_t nCol, size_t nSec, float  sigma);
template
void filter_log<double>(double* const array, size_t nRow, size_t nCol, size_t nSec, double sigma);

} // namespace gem

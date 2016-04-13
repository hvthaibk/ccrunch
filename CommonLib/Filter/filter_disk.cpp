/***********************************************************************
 *  File:       filter_disk.cpp
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
 * Disk
 *************************/

// 1D
template <typename T>
void filter_disk(T* const array, size_t nRow, T radius, bool norm)
{
    assert(array != NULL);
    assert(nRow > 0);
    assert(radius > 0);

    T        cRow = ((T) nRow - 1) / 2;
    T        distToCenter;
    size_t   iRow;

#ifndef NDEBUG
    if (radius > cRow) {
        std::cout << "filter_disk(): "
                  << "radius exceeds filter's boundary"
                  << std::endl;
    }
#endif

    #pragma omp parallel for private(iRow,distToCenter)
    for (iRow = 0; iRow < nRow; iRow++) {
        distToCenter = std::abs((T) iRow - cRow);

        array[iRow] = (T) ((distToCenter <= radius) ? 1 : 0);
    }

    if (norm) {
        T    arraySum = array_reduce_sum(array, nRow);
        array_math_div(array, arraySum, nRow);
    }
}

// instantiation
template
void filter_disk<float >(float*  const array, size_t nRow, float  radius, bool norm);
template
void filter_disk<double>(double* const array, size_t nRow, double radius, bool norm);

// 2D
template <typename T>
void filter_disk(T* const array, size_t nRow, size_t nCol, T radius, bool norm)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0);
    assert(radius > 0);

    T        cRow = ((T) nRow - 1) / 2;
    T        cCol = ((T) nCol - 1) / 2;
    T        dist1, distToCenter;
    size_t   indx1;
    size_t   iRow, iCol;

#ifndef NDEBUG
    if ((radius > cRow) || (radius > cCol)) {
        std::cout << "filter_disk(): "
                  << "radius exceeds filter's boundary"
                  << std::endl;
    }
#endif

    #pragma omp parallel for private(iRow,dist1,indx1,iCol,distToCenter)
    for (iRow = 0; iRow < nRow; iRow++) {
        dist1 = distSqr((T)iRow,cRow);
        indx1 = iRow*nCol;

        for (iCol = 0; iCol < nCol; iCol++) {
            distToCenter = std::sqrt(dist1+distSqr((T)iCol,cCol));

            array[indx1+iCol] = (T) ((distToCenter <= radius) ? 1 : 0);
        }
    }

    if (norm) {
        T    arraySum = array_reduce_sum(array, nRow*nCol);
        array_math_div(array, arraySum, nRow*nCol);
    }
}

// instantiation
template
void filter_disk<float >(float*  const array, size_t nRow, size_t nCol, float  radius, bool norm);
template
void filter_disk<double>(double* const array, size_t nRow, size_t nCol, double radius, bool norm);

// 3D
template <typename T>
void filter_disk(T* const array, size_t nRow, size_t nCol, size_t nSec, T radius, bool norm)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);
    assert(radius > 0);

    T        cRow = ((T) nRow - 1) / 2;
    T        cCol = ((T) nCol - 1) / 2;
    T        cSec = ((T) nSec - 1) / 2;
    T        dist1, dist2, distToCenter;
    size_t   indx1, indx2;
    size_t   iRow, iCol, iSec;

#ifndef NDEBUG
    if ((radius > cRow) || (radius > cCol) || (radius > cSec)) {
        std::cout << "filter_disk(): "
                  << "radius exceeds filter's boundary"
                  << std::endl;
    }
#endif

    #pragma omp parallel for private(iRow,dist1,indx1,iCol,dist2,indx2,iSec,distToCenter)
    for (iRow = 0; iRow < nRow; iRow++) {
        dist1 = distSqr((T)iRow,cRow);
        indx1 = iRow*nCol*nSec;

        for (iCol = 0; iCol < nCol; iCol++) {
            dist2 = dist1 + distSqr((T)iCol,cCol);
            indx2 = indx1 + iCol*nSec;

            for (iSec = 0; iSec < nSec; iSec++) {
                distToCenter = std::sqrt(dist2+distSqr((T)iSec,cSec));

                array[indx2+iSec] = (T) ((distToCenter <= radius) ? 1 : 0);
            }
        }
    }

    if (norm) {
        T    arraySum = array_reduce_sum(array, nRow*nCol*nSec);
        array_math_div(array, arraySum, nRow*nCol*nSec);
    }
}

// instantiation
template
void filter_disk<float >(float*  const array, size_t nRow, size_t nCol, size_t nSec, float  radius, bool norm);
template
void filter_disk<double>(double* const array, size_t nRow, size_t nCol, size_t nSec, double radius, bool norm);

} // namespace gem

/***********************************************************************
 *  File:       array_linearalgebra.cpp
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
 * Linear algebra
 ****************************************/

// matrix multiplication
template <typename T>
void array_linalg_matmul(const T* const arraySrc1, size_t nRowSrc1, size_t nColSrc1,
                         const T* const arraySrc2, size_t nRowSrc2, size_t nColSrc2,
                               T* const arrayDst)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc1 > 0);
    assert(nColSrc1 > 0);
    assert(nRowSrc2 > 0);
    assert(nColSrc2 > 0);

    require(nColSrc1 == nRowSrc2, "array_linalg_matmul: incompatible size");
    for (size_t iRow = 0, i = 0; iRow < nRowSrc1; iRow++) {
        for (size_t iCol = 0; iCol < nColSrc2; iCol++, i++) {
            arrayDst[i] = 0;
            for (size_t iTmp = 0; iTmp < nColSrc1; iTmp++) {
                arrayDst[i] += arraySrc1[sub2ind(iRow,iTmp,nColSrc1)] *
                               arraySrc2[sub2ind(iTmp,iCol,nColSrc2)];
            }
        }
    }
}

// instantiation
template
void array_linalg_matmul<float >(const float*  const arraySrc1, size_t nRowSrc1, size_t nColSrc1,
                                 const float*  const arraySrc2, size_t nRowSrc2, size_t nColSrc2,
                                       float*  const arrayDst);
template
void array_linalg_matmul<double>(const double* const arraySrc1, size_t nRowSrc1, size_t nColSrc1,
                                 const double* const arraySrc2, size_t nRowSrc2, size_t nColSrc2,
                                       double* const arrayDst);

// product
template <typename T>
void array_linalg_prodvec2mat(T* const arrayDst, const T* const arrayRow, size_t nRow,
                                                 const T* const arrayCol, size_t nCol)
{
    assert(arrayDst != NULL);
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(nRow > 0);
    assert(nCol > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        for (size_t iCol = 0; iCol < nCol; iCol++) {
                arrayDst[iRow*nCol+iCol] = arrayRow[iRow]*arrayCol[iCol];
        }
    }
}

template <typename T>
void array_linalg_prodvec2mat(T* const arrayDst, const T* const arrayRow, size_t nRow,
                                                 const T* const arrayCol, size_t nCol,
                                                 const T* const arraySec, size_t nSec)
{
    assert(arrayDst != NULL);
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(arraySec != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
    assert(nSec > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        for (size_t iCol = 0; iCol < nCol; iCol++) {
            for (size_t iSec = 0; iSec < nSec; iSec++) {
                arrayDst[iRow*nCol*nSec+iCol*nSec+iSec] =
                    arrayRow[iRow]*arrayCol[iCol]*arraySec[iSec];
            }
        }
    }
}

// instantiation
template
void array_linalg_prodvec2mat<float >(float*  const arrayDst, const float*  const arrayRow, size_t nRow,
                                                              const float*  const arrayCol, size_t nCol);
template
void array_linalg_prodvec2mat<double>(double* const arrayDst, const double* const arrayRow, size_t nRow,
                                                              const double* const arrayCol, size_t nCol);
template
void array_linalg_prodvec2mat<float >(float*  const arrayDst, const float*  const arrayRow, size_t nRow,
                                                              const float*  const arrayCol, size_t nCol,
                                                              const float*  const arraySec, size_t nSec);
template
void array_linalg_prodvec2mat<double>(double* const arrayDst, const double* const arrayRow, size_t nRow,
                                                              const double* const arrayCol, size_t nCol,
                                                              const double* const arraySec, size_t nSec);

// blas axpy
template <typename T>
void array_linalg_axpy(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] += value*arraySrc[iRow];
    }
}

// instantiation
template
void array_linalg_axpy<float >(float*  const arrayDst, const float*  const arraySrc, float  value, size_t nRow);
template
void array_linalg_axpy<double>(double* const arrayDst, const double* const arraySrc, double value, size_t nRow);

} // namespace gem

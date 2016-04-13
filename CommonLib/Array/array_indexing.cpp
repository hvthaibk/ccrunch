/***********************************************************************
 *  File:       array_indexing.cpp
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
 * Indexing
 ****************************************/

// 1D index -> 2D sub-indices
#ifndef NDEBUG
void array_index_ind2sub(const size_t* const arrayInd,
                               size_t* const arrayRow,
                               size_t* const arrayCol,
                               size_t        length,
                               size_t        nRow, size_t nCol)
#else
void array_index_ind2sub(const size_t* const arrayInd,
                               size_t* const arrayRow,
                               size_t* const arrayCol,
                               size_t        length,
                               size_t            , size_t nCol)
#endif
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(arrayInd != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0);

    #pragma omp parallel for
    for (size_t i = 0; i < length; i++) {
        assert(arrayInd[i] < nRow*nCol);

        arrayRow[i] = arrayInd[i] / nCol;
        arrayCol[i] = arrayInd[i] % nCol;
    }
}

// 1D index -> 3D sub-indices
#ifndef NDEBUG
void array_index_ind2sub(const size_t* const arrayInd,
                               size_t* const arrayRow,
                               size_t* const arrayCol,
                               size_t* const arraySec,
                               size_t        length,
                               size_t        nRow, size_t nCol, size_t nSec)
#else
void array_index_ind2sub(const size_t* const arrayInd,
                               size_t* const arrayRow,
                               size_t* const arrayCol,
                               size_t* const arraySec,
                               size_t        length,
                               size_t            , size_t nCol, size_t nSec)
#endif
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(arraySec != NULL);
    assert(arrayInd != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    size_t    nColSec = nCol*nSec;

    #pragma omp parallel for
    for (size_t i = 0; i < length; i++) {
        assert(arrayInd[i] < nRow*nCol*nSec);

        arrayRow[i] = arrayInd[i] / nColSec;
        arrayCol[i] = (arrayInd[i] - arrayRow[i]*nColSec) / nSec;
        arraySec[i] = arrayInd[i] % nSec;
    }
}

// 2D sub-indices -> 1D index
#ifndef NDEBUG
void array_index_sub2ind(const size_t* const arrayRow,
                         const size_t* const arrayCol,
                               size_t* const arrayInd,
                               size_t        length,
                               size_t        nRow, size_t nCol)
#else
void array_index_sub2ind(const size_t* const arrayRow,
                         const size_t* const arrayCol,
                               size_t* const arrayInd,
                               size_t        length,
                               size_t            , size_t nCol)
#endif
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(arrayInd != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0);

    #pragma omp parallel for
    for (size_t i = 0; i < length; i++) {
        assert(arrayRow[i] < nRow);
        assert(arrayCol[i] < nCol);

        arrayInd[i] = arrayRow[i]*nCol + arrayCol[i];
    }
}

// 3D sub-indices -> 1D index
#ifndef NDEBUG
void array_index_sub2ind(const size_t* const arrayRow,
                         const size_t* const arrayCol,
                         const size_t* const arraySec,
                               size_t* const arrayInd,
                               size_t        length,
                               size_t        nRow, size_t nCol, size_t nSec)
#else
void array_index_sub2ind(const size_t* const arrayRow,
                         const size_t* const arrayCol,
                         const size_t* const arraySec,
                               size_t* const arrayInd,
                               size_t        length,
                               size_t            , size_t nCol, size_t nSec)
#endif
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(arraySec != NULL);
    assert(arrayInd != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    #pragma omp parallel for
    for (size_t i = 0; i < length; i++) {
        assert(arrayRow[i] < nRow);
        assert(arrayCol[i] < nCol);
        assert(arraySec[i] < nSec);

        arrayInd[i] = (arrayRow[i]*nCol + arrayCol[i])*nSec + arraySec[i];
    }
}

// column-major -> row-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow)
template <typename T1, typename T2>
void array_index_col2row(const T1* const arrayCol, T2* const arrayRow,
                         size_t nRow, size_t nCol)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(nRow > 0 && nCol > 0);

    size_t    iCol, iRow;

#ifdef __GEM_USE_OPENMP__
    size_t    indxColTmp;

    #pragma omp parallel for private(iCol,indxColTmp,iRow)
    for (iCol = 0; iCol < nCol; iCol++) {
        indxColTmp = iCol*nRow;

        for (iRow = 0; iRow < nRow; iRow++) {
            arrayRow[iRow*nCol+iCol] = (T2) arrayCol[indxColTmp+iRow];
        }
    }
#else
    size_t    indxCol;

    for (iCol = 0, indxCol = 0; iCol < nCol; iCol++) {
        for (iRow = 0; iRow < nRow; iRow++, indxCol++) {
            arrayRow[iRow*nCol+iCol] = (T2) arrayCol[indxCol];
        }
    }
#endif
}

// instantiation
template
void array_index_col2row<int32_t ,int32_t >(const int32_t*  const arrayCol, int32_t*  const arrayRow,
                                            size_t nRow, size_t nCol);
template
void array_index_col2row<uint32_t,uint32_t>(const uint32_t* const arrayCol, uint32_t* const arrayRow,
                                            size_t nRow, size_t nCol);
template
void array_index_col2row<int64_t ,int64_t >(const int64_t*  const arrayCol, int64_t*  const arrayRow,
                                            size_t nRow, size_t nCol);
template
void array_index_col2row<uint64_t,uint64_t>(const uint64_t* const arrayCol, uint64_t* const arrayRow,
                                            size_t nRow, size_t nCol);
template
void array_index_col2row<float   ,float   >(const float*    const arrayCol, float*    const arrayRow,
                                            size_t nRow, size_t nCol);
template
void array_index_col2row<float   ,double  >(const float*    const arrayCol, double*   const arrayRow,
                                            size_t nRow, size_t nCol);
template
void array_index_col2row<double  ,float   >(const double*   const arrayCol, float*    const arrayRow,
                                            size_t nRow, size_t nCol);
template
void array_index_col2row<double  ,double  >(const double*   const arrayCol, double*   const arrayRow,
                                            size_t nRow, size_t nCol);

// column-major -> row-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow), nSec = dim3(arrayRow)
template <typename T1, typename T2>
void array_index_col2row(const T1* const arrayCol, T2* const arrayRow,
                         size_t nRow, size_t nCol, size_t nSec)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    size_t    iSec, iCol, iRow;

#ifdef __GEM_USE_OPENMP__
    size_t    nRowCol = nRow*nCol;
    size_t    nColSec = nCol*nSec;
    size_t    indxRowTmp, indxColTmp1, indxColTmp2;

    #pragma omp parallel for private(iSec,indxColTmp1,iCol,indxRowTmp,indxColTmp2,iRow)
    for (iSec = 0; iSec < nSec; iSec++) {
        indxColTmp1 = iSec*nRowCol;

        for (iCol = 0; iCol < nCol; iCol++) {
            indxRowTmp  = iCol*nSec + iSec;
            indxColTmp2 = indxColTmp1 + iCol*nRow;

            for (iRow = 0; iRow < nRow; iRow++) {
                arrayRow[iRow*nColSec+indxRowTmp] = (T2) arrayCol[indxColTmp2+iRow];
            }
        }
    }
#else
    size_t    indxCol, indxRowTmp;
    size_t    nColSec = nCol*nSec;

    for (iSec = 0, indxCol = 0; iSec < nSec; iSec++) {
        for (iCol = 0; iCol < nCol; iCol++) {
            indxRowTmp = iCol*nSec + iSec;

            for (iRow = 0; iRow < nRow; iRow++, indxCol++) {
                arrayRow[iRow*nColSec+indxRowTmp] = (T2) arrayCol[indxCol];
            }
        }
    }
#endif
}

// instantiation
template
void array_index_col2row<int32_t ,int32_t >(const int32_t*  const arrayCol, int32_t*  const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<uint32_t,uint32_t>(const uint32_t* const arrayCol, uint32_t* const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<int64_t ,int64_t >(const int64_t*  const arrayCol, int64_t*  const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<uint64_t,uint64_t>(const uint64_t* const arrayCol, uint64_t* const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<int8_t  ,float   >(const int8_t*   const arrayCol, float*    const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<int8_t  ,double  >(const int8_t*   const arrayCol, double*   const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<uint16_t,float   >(const uint16_t* const arrayCol, float*    const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<uint16_t,double  >(const uint16_t* const arrayCol, double*   const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<int16_t ,float   >(const int16_t*  const arrayCol, float*    const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<int16_t ,double  >(const int16_t*  const arrayCol, double*   const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<float   ,float   >(const float*    const arrayCol, float*    const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<float   ,double  >(const float*    const arrayCol, double*   const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<double  ,float   >(const double*   const arrayCol, float*    const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_col2row<double  ,double  >(const double*   const arrayCol, double*   const arrayRow,
                                            size_t nRow, size_t nCol, size_t nSec);


// row-major -> column-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow)
template <typename T1, typename T2>
void array_index_row2col(const T1* const arrayRow, T2* const arrayCol,
                         size_t nRow, size_t nCol)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(nRow > 0 && nCol > 0);

    size_t    iCol, iRow;

#ifdef __GEM_USE_OPENMP__
    size_t    indxColTmp;

    #pragma omp parallel for private(iCol,indxColTmp,iRow)
    for (iCol = 0; iCol < nCol; iCol++) {
        indxColTmp = iCol*nRow;

        for (iRow = 0; iRow < nRow; iRow++) {
            arrayCol[indxColTmp+iRow] = (T2) arrayRow[iRow*nCol+iCol];
        }
    }
#else
    size_t    indxCol;

    for (iCol = 0, indxCol = 0; iCol < nCol; iCol++) {
        for (iRow = 0; iRow < nRow; iRow++, indxCol++) {
            arrayCol[indxCol] = (T2) arrayRow[iRow*nCol+iCol];
        }
    }
#endif
}

// instantiation
template
void array_index_row2col<int32_t ,int32_t >(const int32_t*  const arrayRow, int32_t*  const arrayCol,
                                            size_t nRow, size_t nCol);
template
void array_index_row2col<uint32_t,uint32_t>(const uint32_t* const arrayRow, uint32_t* const arrayCol,
                                            size_t nRow, size_t nCol);
template
void array_index_row2col<int64_t ,int64_t >(const int64_t*  const arrayRow, int64_t*  const arrayCol,
                                            size_t nRow, size_t nCol);
template
void array_index_row2col<uint64_t,uint64_t>(const uint64_t* const arrayRow, uint64_t* const arrayCol,
                                            size_t nRow, size_t nCol);
template
void array_index_row2col<float   ,float   >(const float*    const arrayRow, float*    const arrayCol,
                                            size_t nRow, size_t nCol);
template
void array_index_row2col<float   ,double  >(const float*    const arrayRow, double*   const arrayCol,
                                            size_t nRow, size_t nCol);
template
void array_index_row2col<double  ,float   >(const double*   const arrayRow, float*    const arrayCol,
                                            size_t nRow, size_t nCol);
template
void array_index_row2col<double  ,double  >(const double*   const arrayRow, double*   const arrayCol,
                                            size_t nRow, size_t nCol);

// row-major -> column-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow), nSec = dim3(arrayRow)
template <typename T1, typename T2>
void array_index_row2col(const T1* const arrayRow, T2* const arrayCol,
                         size_t nRow, size_t nCol, size_t nSec)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    size_t    iSec, iCol, iRow;

#ifdef __GEM_USE_OPENMP__
    size_t    nRowCol = nRow*nCol;
    size_t    nColSec = nCol*nSec;
    size_t    indxRowTmp, indxColTmp1, indxColTmp2;

    #pragma omp parallel for private(iSec,indxColTmp1,iCol,indxRowTmp,indxColTmp2,iRow)
    for (iSec = 0; iSec < nSec; iSec++) {
        indxColTmp1 = iSec*nRowCol;

        for (iCol = 0; iCol < nCol; iCol++) {
            indxRowTmp  = iCol*nSec + iSec;
            indxColTmp2 = indxColTmp1 + iCol*nRow;

            for (iRow = 0; iRow < nRow; iRow++) {
                arrayCol[indxColTmp2+iRow] = (T2) arrayRow[iRow*nColSec+indxRowTmp];
            }
        }
    }
#else
    size_t    indxCol, indxRowTmp;
    size_t    nColSec = nCol*nSec;

    for (iSec = 0, indxCol = 0; iSec < nSec; iSec++) {
        for (iCol = 0; iCol < nCol; iCol++) {
            indxRowTmp = iCol*nSec + iSec;

            for (iRow = 0; iRow < nRow; iRow++, indxCol++) {
                arrayCol[indxCol] = (T2) arrayRow[iRow*nColSec+indxRowTmp];
            }
        }
    }
#endif
}

// instantiation
template
void array_index_row2col<int32_t ,int32_t >(const int32_t*  const arrayRow, int32_t*  const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<uint32_t,uint32_t>(const uint32_t* const arrayRow, uint32_t* const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<int64_t ,int64_t >(const int64_t*  const arrayRow, int64_t*  const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<uint64_t,uint64_t>(const uint64_t* const arrayRow, uint64_t* const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<float   ,int8_t  >(const float*    const arrayRow, int8_t*   const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<double  ,int8_t  >(const double*   const arrayRow, int8_t*   const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<float   ,uint16_t>(const float*    const arrayRow, uint16_t* const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<double  ,uint16_t>(const double*   const arrayRow, uint16_t* const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<float   ,int16_t >(const float*    const arrayRow, int16_t*  const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<double  ,int16_t >(const double*   const arrayRow, int16_t*  const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<float   ,float   >(const float*    const arrayRow, float*    const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<float   ,double  >(const float*    const arrayRow, double*   const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<double  ,float   >(const double*   const arrayRow, float*    const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);
template
void array_index_row2col<double  ,double  >(const double*   const arrayRow, double*   const arrayCol,
                                            size_t nRow, size_t nCol, size_t nSec);

} // namespace gem

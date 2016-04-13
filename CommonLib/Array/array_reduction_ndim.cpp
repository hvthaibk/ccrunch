/***********************************************************************
 *  File:       array_reduction_ndim.cpp
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
 * Reduction NDim
 ****************************************/

template <typename T>
void array_reducendim_max(const T* const arraySrc, size_t nRow, size_t nCol,
                                T* const arrayDst,
                          eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);

    size_t    iRow, iCol;
    size_t    iSrc;

    switch (dim) {
        case REDUCE_NDIM_ROW:
            for (iCol = 0; iCol < nCol; iCol++) {
                arrayDst[iCol] = arraySrc[iCol];
            }

            for (iRow = 1, iSrc = nCol; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++, iSrc++) {
                    arrayDst[iCol] = max(arraySrc[iSrc], arrayDst[iCol]);
                }
            }
            break;
        case REDUCE_NDIM_COL:
            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                arrayDst[iRow] = arraySrc[iSrc++];

                for (iCol = 1; iCol < nCol; iCol++, iSrc++) {
                    arrayDst[iRow] = max(arraySrc[iSrc], arrayDst[iRow]);
                }
            }
            break;
        default:
            ERROR("array_reducendim_max", "unsupported dimension mode");
    }
}

// instantiation
template
void array_reducendim_max<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol,
                                          int32_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_max<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol,
                                          uint32_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_max<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol,
                                          int64_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_max<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol,
                                          uint64_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_max<float   >(const float*    const arraySrc, size_t nRow, size_t nCol,
                                          float*    const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_max<double  >(const double*   const arraySrc, size_t nRow, size_t nCol,
                                          double*   const arrayDst,
                                    eReduceNDim dim);

template <typename T>
void array_reducendim_min(const T* const arraySrc, size_t nRow, size_t nCol,
                                T* const arrayDst,
                          eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);

    size_t    iRow, iCol;
    size_t    iSrc;

    switch (dim) {
        case REDUCE_NDIM_ROW:
            for (iCol = 0; iCol < nCol; iCol++) {
                arrayDst[iCol] = arraySrc[iCol];
            }

            for (iRow = 1, iSrc = nCol; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++, iSrc++) {
                    arrayDst[iCol] = min(arraySrc[iSrc], arrayDst[iCol]);
                }
            }
            break;
        case REDUCE_NDIM_COL:
            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                arrayDst[iRow] = arraySrc[iSrc++];

                for (iCol = 1; iCol < nCol; iCol++, iSrc++) {
                    arrayDst[iRow] = min(arraySrc[iSrc], arrayDst[iRow]);
                }
            }
            break;
        default:
            ERROR("array_reducendim_min", "unsupported dimension mode");
    }
}

// instantiation
template
void array_reducendim_min<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol,
                                          int32_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_min<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol,
                                          uint32_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_min<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol,
                                          int64_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_min<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol,
                                          uint64_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_min<float   >(const float*    const arraySrc, size_t nRow, size_t nCol,
                                          float*    const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_min<double  >(const double*   const arraySrc, size_t nRow, size_t nCol,
                                          double*   const arrayDst,
                                    eReduceNDim dim);

template <typename T>
void array_reducendim_sum(const T* const arraySrc, size_t nRow, size_t nCol,
                                T* const arrayDst,
                          eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);

    size_t    iRow, iCol;
    size_t    iSrc;

    switch (dim) {
        case REDUCE_NDIM_ROW:
            array_memset(arrayDst, 0, nCol);

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++, iSrc++) {
                    arrayDst[iCol] += arraySrc[iSrc];
                }
            }
            break;
        case REDUCE_NDIM_COL:
            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                arrayDst[iRow] = 0;

                for (iCol = 0; iCol < nCol; iCol++, iSrc++) {
                    arrayDst[iRow] += arraySrc[iSrc];
                }
            }
            break;
        default:
            ERROR("array_reducendim_sum", "unsupported dimension mode");
    }
}

// instantiation
template
void array_reducendim_sum<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol,
                                          int32_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_sum<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol,
                                          uint32_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_sum<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol,
                                          int64_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_sum<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol,
                                          uint64_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_sum<float   >(const float*    const arraySrc, size_t nRow, size_t nCol,
                                          float*    const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_sum<double  >(const double*   const arraySrc, size_t nRow, size_t nCol,
                                          double*   const arrayDst,
                                    eReduceNDim dim);

template <typename T>
void array_reducendim_sum2(const T* const arraySrc, size_t nRow, size_t nCol,
                                 T* const arrayDst,
                           eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);

    size_t    iRow, iCol;
    size_t    iSrc;

    switch (dim) {
        case REDUCE_NDIM_ROW:
            array_memset(arrayDst, 0, nCol);

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++, iSrc++) {
                    arrayDst[iCol] += pow2(arraySrc[iSrc]);
                }
            }
            break;
        case REDUCE_NDIM_COL:
            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                arrayDst[iRow] = 0;

                for (iCol = 0; iCol < nCol; iCol++, iSrc++) {
                    arrayDst[iRow] += pow2(arraySrc[iSrc]);
                }
            }
            break;
        default:
            ERROR("array_reducendim_sum2", "unsupported dimension mode");
    }
}

// instantiation
template
void array_reducendim_sum2<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol,
                                           int32_t*  const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_sum2<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol,
                                           uint32_t* const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_sum2<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol,
                                           int64_t*  const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_sum2<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol,
                                           uint64_t* const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_sum2<float   >(const float*    const arraySrc, size_t nRow, size_t nCol,
                                           float*    const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_sum2<double  >(const double*   const arraySrc, size_t nRow, size_t nCol,
                                           double*   const arrayDst,
                                     eReduceNDim dim);

template <typename T>
void array_reducendim_mean(const T* const arraySrc, size_t nRow, size_t nCol,
                                 T* const arrayDst,
                           eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);

    array_reducendim_sum(arraySrc, nRow, nCol,
                         arrayDst,
                         dim);

    switch (dim) {
        case REDUCE_NDIM_ROW:
            array_math_div(arrayDst, (T) nRow, nCol);
            break;
        case REDUCE_NDIM_COL:
            array_math_div(arrayDst, (T) nCol, nRow);
            break;
        default:
            ERROR("array_reducendim_mean", "unsupported dimension mode");
    }
}

// instantiation
template
void array_reducendim_mean<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol,
                                           int32_t*  const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_mean<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol,
                                           uint32_t* const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_mean<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol,
                                           int64_t*  const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_mean<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol,
                                           uint64_t* const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_mean<float   >(const float*    const arraySrc, size_t nRow, size_t nCol,
                                           float*    const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_mean<double  >(const double*   const arraySrc, size_t nRow, size_t nCol,
                                           double*   const arrayDst,
                                     eReduceNDim dim);

template <typename T>
void array_reducendim_std (const T* const arraySrc, size_t nRow, size_t nCol,
                                 T* const arrayDst,
                           eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);

    T*        arrayMean = NULL;
    size_t    iRow, iCol;
    size_t    iSrc;
    T         sum;

    switch (dim) {
        case REDUCE_NDIM_ROW:
            array_new(arrayMean, nCol);
            break;
        case REDUCE_NDIM_COL:
            array_new(arrayMean, nRow);
            break;
        default:
            ERROR("array_reducendim_std", "unsupported dimension mode");
    }

    array_reducendim_mean(arraySrc, nRow, nCol,
                          arrayMean,
                          dim);

    switch (dim) {
        case REDUCE_NDIM_ROW:
            array_memset(arrayDst, 0, nCol);

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++, iSrc++) {
                    arrayDst[iCol] += pow2(arraySrc[iSrc]-arrayMean[iCol]);
                }
            }

            array_math_div (arrayDst, (T) nRow, nCol);
            array_math_sqrt(arrayDst, nCol);
            break;
        case REDUCE_NDIM_COL:
            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                arrayDst[iRow] = 0;

                for (iCol = 0; iCol < nCol; iCol++) {
                    sum += pow2(arraySrc[iSrc]-arrayMean[iRow]);
                }
            }

            array_math_div (arrayDst, (T) nCol, nRow);
            array_math_sqrt(arrayDst, nRow);
            break;
        default:
            ERROR("array_reducendim_std", "unsupported dimension mode");
    }

    array_delete(arrayMean);
}

// instantiation
template
void array_reducendim_std <int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol,
                                           int32_t*  const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_std <uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol,
                                           uint32_t* const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_std <int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol,
                                           int64_t*  const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_std <uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol,
                                           uint64_t* const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_std <float   >(const float*    const arraySrc, size_t nRow, size_t nCol,
                                           float*    const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_std <double  >(const double*   const arraySrc, size_t nRow, size_t nCol,
                                           double*   const arrayDst,
                                     eReduceNDim dim);

template <typename T>
void array_reducendim_max(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                T* const arrayDst,
                          eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
    assert(nSec > 0);

    size_t    iRow, iCol, iSec;
    size_t    iSrc, iDst;

    switch (dim) {
        case REDUCE_NDIM_ROW:
            for (iCol = 0, iSrc = 0, iDst = 0; iCol < nCol; iCol++) {
                for (iSec = 0; iSec < nSec; iSec++, iSrc++, iDst++) {
                    arrayDst[iDst] = arraySrc[iSrc];
                }
            }

            for (iRow = 1, iSrc = nCol*nSec; iRow < nRow; iRow++) {
                for (iCol = 0, iDst = 0; iCol < nCol; iCol++) {
                    for (iSec = 0; iSec < nSec; iSec++, iSrc++, iDst++) {
                        arrayDst[iDst] = max(arraySrc[iSrc], arrayDst[iDst]);
                    }
                }
            }
            break;
        case REDUCE_NDIM_COL:
            for (iRow = 0, iDst = 0; iRow < nRow; iRow++) {
                for (iSec = 0; iSec < nSec; iSec++, iDst++) {
                    arrayDst[iDst] = arraySrc[sub2ind(iRow,0,iSec,nCol,nSec)];
                }
            }

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        iDst = sub2ind(iRow,iSec,nSec);
                        arrayDst[iDst] = max(arraySrc[iSrc], arrayDst[iDst]);
                    }
                }
            }
            break;
        case REDUCE_NDIM_SEC:
            for (iRow = 0, iSrc = 0, iDst = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++, iDst++) {
                    arrayDst[iDst] = arraySrc[iSrc++];

                    for (iSec = 1; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[iDst] = max(arraySrc[iSrc], arrayDst[iDst]);
                    }
                }
            }
            break;
        default:
            ERROR("array_reducendim_max", "unsupported dimension mode");
    }
}

// instantiation
template
void array_reducendim_max<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          int32_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_max<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          uint32_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_max<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          int64_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_max<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          uint64_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_max<float   >(const float*    const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          float*    const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_max<double  >(const double*   const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          double*   const arrayDst,
                                    eReduceNDim dim);

template <typename T>
void array_reducendim_min(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                T* const arrayDst,
                          eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
    assert(nSec > 0);

    size_t    iRow, iCol, iSec;
    size_t    iSrc, iDst;

    switch (dim) {
        case REDUCE_NDIM_ROW:
            for (iCol = 0, iSrc = 0, iDst = 0; iCol < nCol; iCol++) {
                for (iSec = 0; iSec < nSec; iSec++, iSrc++, iDst++) {
                    arrayDst[iDst] = arraySrc[iSrc];
                }
            }

            for (iRow = 1, iSrc = nCol*nSec; iRow < nRow; iRow++) {
                for (iCol = 0, iDst = 0; iCol < nCol; iCol++) {
                    for (iSec = 0; iSec < nSec; iSec++, iSrc++, iDst++) {
                        arrayDst[iDst] = min(arraySrc[iSrc], arrayDst[iDst]);
                    }
                }
            }
            break;
        case REDUCE_NDIM_COL:
            for (iRow = 0, iDst = 0; iRow < nRow; iRow++) {
                for (iSec = 0; iSec < nSec; iSec++, iDst++) {
                    arrayDst[iDst] = arraySrc[sub2ind(iRow,0,iSec,nCol,nSec)];
                }
            }

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        iDst = sub2ind(iRow,iSec,nSec);
                        arrayDst[iDst] = min(arraySrc[iSrc], arrayDst[iDst]);
                    }
                }
            }
            break;
        case REDUCE_NDIM_SEC:
            for (iRow = 0, iSrc = 0, iDst = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++, iDst++) {
                    arrayDst[iDst] = arraySrc[iSrc++];

                    for (iSec = 1; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[iDst] = min(arraySrc[iSrc], arrayDst[iDst]);
                    }
                }
            }
            break;
        default:
            ERROR("array_reducendim_min", "unsupported dimension mode");
    }
}

// instantiation
template
void array_reducendim_min<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          int32_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_min<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          uint32_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_min<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          int64_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_min<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          uint64_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_min<float   >(const float*    const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          float*    const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_min<double  >(const double*   const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          double*   const arrayDst,
                                    eReduceNDim dim);

template <typename T>
void array_reducendim_sum(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                T* const arrayDst,
                          eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
    assert(nSec > 0);

    size_t    iRow, iCol, iSec;
    size_t    iSrc, iDst;

    switch (dim) {
        case REDUCE_NDIM_ROW:
            array_memset(arrayDst, 0, nCol*nSec);

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[sub2ind(iCol,iSec,nSec)] += arraySrc[iSrc];
                    }
                }
            }
            break;
        case REDUCE_NDIM_COL:
            array_memset(arrayDst, 0, nRow*nSec);

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[sub2ind(iRow,iSec,nSec)] += arraySrc[iSrc];
                    }
                }
            }
            break;
        case REDUCE_NDIM_SEC:
            for (iRow = 0, iDst = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++, iDst++) {
                    arrayDst[iDst] = 0;

                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[iDst] += arraySrc[iSrc];
                    }
                }
            }
            break;
        default:
            ERROR("array_reducendim_sum", "unsupported dimension mode");
    }
}

// instantiation
template
void array_reducendim_sum<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          int32_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_sum<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          uint32_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_sum<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          int64_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_sum<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          uint64_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_sum<float   >(const float*    const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          float*    const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_sum<double  >(const double*   const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          double*   const arrayDst,
                                    eReduceNDim dim);

template <typename T>
void array_reducendim_sum2(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                 T* const arrayDst,
                           eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
    assert(nSec > 0);

    size_t    iRow, iCol, iSec;
    size_t    iSrc, iDst;

    switch (dim) {
        case REDUCE_NDIM_ROW:
            array_memset(arrayDst, 0, nCol*nSec);

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[sub2ind(iCol,iSec,nSec)] += pow2(arraySrc[iSrc]);
                    }
                }
            }
            break;
        case REDUCE_NDIM_COL:
            array_memset(arrayDst, 0, nRow*nSec);

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[sub2ind(iRow,iSec,nSec)] += pow2(arraySrc[iSrc]);
                    }
                }
            }
            break;
        case REDUCE_NDIM_SEC:
            for (iRow = 0, iDst = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++, iDst++) {
                    arrayDst[iDst] = 0;

                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[iDst] += pow2(arraySrc[iSrc]);
                    }
                }
            }
            break;
        default:
            ERROR("array_reducendim_sum2", "unsupported dimension mode");
    }
}

// instantiation
template
void array_reducendim_sum2<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           int32_t*  const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_sum2<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           uint32_t* const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_sum2<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           int64_t*  const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_sum2<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           uint64_t* const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_sum2<float   >(const float*    const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           float*    const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_sum2<double  >(const double*   const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           double*   const arrayDst,
                                     eReduceNDim dim);

template <typename T>
void array_reducendim_mean(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                 T* const arrayDst,
                           eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
    assert(nSec > 0);

    array_reducendim_sum(arraySrc, nRow, nCol, nSec,
                         arrayDst,
                         dim);

    switch (dim) {
        case REDUCE_NDIM_ROW:
            array_math_div(arrayDst, (T) nRow, nCol*nSec);
            break;
        case REDUCE_NDIM_COL:
            array_math_div(arrayDst, (T) nCol, nRow*nSec);
            break;
        case REDUCE_NDIM_SEC:
            array_math_div(arrayDst, (T) nSec, nRow*nCol);
            break;
        default:
            ERROR("array_reducendim_mean", "unsupported dimension mode");
    }
}

// instantiation
template
void array_reducendim_mean<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           int32_t*  const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_mean<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           uint32_t* const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_mean<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           int64_t*  const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_mean<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           uint64_t* const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_mean<float   >(const float*    const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           float*    const arrayDst,
                                     eReduceNDim dim);
template
void array_reducendim_mean<double  >(const double*   const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                           double*   const arrayDst,
                                     eReduceNDim dim);

template <typename T>
void array_reducendim_std(const T* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                T* const arrayDst,
                          eReduceNDim dim)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
    assert(nSec > 0);

    T*        arrayMean = NULL;
    size_t    iRow, iCol, iSec;
    size_t    iSrc, iDst;

    switch (dim) {
        case REDUCE_NDIM_ROW:
            array_new(arrayMean, nCol*nSec);
            break;
        case REDUCE_NDIM_COL:
            array_new(arrayMean, nRow*nSec);
            break;
        case REDUCE_NDIM_SEC:
            array_new(arrayMean, nRow*nCol);
            break;
        default:
            ERROR("array_reducendim_std", "unsupported dimension mode");
    }

    array_reducendim_mean(arraySrc, nRow, nCol, nSec,
                          arrayMean,
                          dim);

    switch (dim) {
        case REDUCE_NDIM_ROW:
            array_memset(arrayDst, 0, nCol*nSec);

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[sub2ind(iCol,iSec,nSec)] += pow2(arraySrc[iSrc]-arrayMean[sub2ind(iCol,iSec,nSec)]);
                    }
                }
            }

            array_math_div (arrayDst, (T) nRow, nCol*nSec);
            array_math_sqrt(arrayDst, nCol*nSec);
            break;
        case REDUCE_NDIM_COL:
            array_memset(arrayDst, 0, nRow*nSec);

            for (iRow = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[sub2ind(iRow,iSec,nSec)] += pow2(arraySrc[iSrc]-arrayMean[sub2ind(iRow,iSec,nSec)]);
                    }
                }
            }

            array_math_div (arrayDst, (T) nCol, nRow*nSec);
            array_math_sqrt(arrayDst, nRow*nSec);
            break;
        case REDUCE_NDIM_SEC:
            for (iRow = 0, iDst = 0, iSrc = 0; iRow < nRow; iRow++) {
                for (iCol = 0; iCol < nCol; iCol++, iDst++) {
                    arrayDst[iDst] = 0;

                    for (iSec = 0; iSec < nSec; iSec++, iSrc++) {
                        arrayDst[iDst] += pow2(arraySrc[iSrc]-arrayMean[iDst]);
                    }
                }
            }

            array_math_div (arrayDst, (T) nSec, nRow*nCol);
            array_math_sqrt(arrayDst, nRow*nCol);
            break;
        default:
            ERROR("array_reducendim_std", "unsupported dimension mode");
    }

    array_delete(arrayMean);
}

// instantiation
template
void array_reducendim_std<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          int32_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_std<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          uint32_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_std<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          int64_t*  const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_std<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          uint64_t* const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_std<float   >(const float*    const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          float*    const arrayDst,
                                    eReduceNDim dim);
template
void array_reducendim_std<double  >(const double*   const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                          double*   const arrayDst,
                                    eReduceNDim dim);

} // namespace gem

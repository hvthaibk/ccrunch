/***********************************************************************
 *  File:       array_permutation.cpp
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
 * Permutation
 ****************************************/

// 2D
template <typename T>
void array_permute(const T* const arraySrc, T* const arrayDst,
                   size_t nRow, size_t nCol)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);

    size_t    iCol, iRow;

#ifdef USE_OPENMP
    size_t    indxTransTmp;

    #pragma omp parallel for private(iCol,indxTransTmp,iRow)
    for (iCol = 0; iCol < nCol; iCol++) {
        indxTransTmp = iCol*nRow;

        for (iRow = 0; iRow < nRow; iRow++) {
            arrayDst[indxTransTmp+iRow] = arraySrc[iRow*nCol+iCol];
        }
    }
#else
    size_t    indxTrans;

    for (iCol = 0, indxTrans = 0; iCol < nCol; iCol++) {
        for (iRow = 0; iRow < nRow; iRow++, indxTrans++) {
            arrayDst[indxTrans] = arraySrc[iRow*nCol+iCol];
        }
    }
#endif
}

// instantiation
template
void array_permute<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                             size_t nRow, size_t nCol);
template
void array_permute<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                             size_t nRow, size_t nCol);
template
void array_permute<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                             size_t nRow, size_t nCol);
template
void array_permute<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                             size_t nRow, size_t nCol);
template
void array_permute<float   >(const float*    const arraySrc, float*    const arrayDst,
                             size_t nRow, size_t nCol);
template
void array_permute<double  >(const double*   const arraySrc, double*   const arrayDst,
                             size_t nRow, size_t nCol);
template
void array_permute<std::complex<float > >(const std::complex<float >* const arraySrc,
                                                std::complex<float >* const arrayDst,
                                                size_t nRow, size_t nCol);
template
void array_permute<std::complex<double> >(const std::complex<double>* const arraySrc,
                                                std::complex<double>* const arrayDst,
                                                size_t nRow, size_t nCol);

// 3D
template <typename T>
void array_permute(const T* const arraySrc, T* const arrayDst,
                   size_t nRow, size_t nCol, size_t nSec,
                   ePermute permute)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    size_t    iSec, iCol, iRow;
    size_t    indxTrans;
    size_t    nColSec = nCol*nSec;

    switch (permute) {
        case PERMUTE3D_123:
            ERROR("array_permute", "no need to permute, memcpy is faster");
            break;
        case PERMUTE3D_132:
            for (iRow = 0, indxTrans = 0; iRow < nRow; iRow++) {
                for (iSec = 0; iSec < nSec; iSec++) {
                    for (iCol = 0; iCol < nCol; iCol++, indxTrans++) {
                        arrayDst[indxTrans] = arraySrc[iRow*nColSec+iCol*nSec+iSec];
                    }
                }
            }
            break;
        case PERMUTE3D_213:
            for (iCol = 0, indxTrans = 0; iCol < nCol; iCol++) {
                for (iRow = 0; iRow < nRow; iRow++) {
                    for (iSec = 0; iSec < nSec; iSec++, indxTrans++) {
                        arrayDst[indxTrans] = arraySrc[iRow*nColSec+iCol*nSec+iSec];
                    }
                }
            }
            break;
        case PERMUTE3D_231:
            for (iCol = 0, indxTrans = 0; iCol < nCol; iCol++) {
                for (iSec = 0; iSec < nSec; iSec++) {
                    for (iRow = 0; iRow < nRow; iRow++, indxTrans++) {
                        arrayDst[indxTrans] = arraySrc[iRow*nColSec+iCol*nSec+iSec];
                    }
                }
            }
            break;
        case PERMUTE3D_312:
            for (iSec = 0, indxTrans = 0; iSec < nSec; iSec++) {
                for (iRow = 0; iRow < nRow; iRow++) {
                    for (iCol = 0; iCol < nCol; iCol++, indxTrans++) {
                        arrayDst[indxTrans] = arraySrc[iRow*nColSec+iCol*nSec+iSec];
                    }
                }
            }
            break;
        case PERMUTE3D_321:
            for (iSec = 0, indxTrans = 0; iSec < nSec; iSec++) {
                for (iCol = 0; iCol < nCol; iCol++) {
                    for (iRow = 0; iRow < nRow; iRow++, indxTrans++) {
                        arrayDst[indxTrans] = arraySrc[iRow*nColSec+iCol*nSec+iSec];
                    }
                }
            }
            break;
        default:
            ERROR("array_permute", "unsupported permutation mode");
    }
}

// instantiation
template
void array_permute<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                             size_t nRow, size_t nCol, size_t nSec,
                             ePermute permute);
template
void array_permute<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                             size_t nRow, size_t nCol, size_t nSec,
                             ePermute permute);
template
void array_permute<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                             size_t nRow, size_t nCol, size_t nSec,
                             ePermute permute);
template
void array_permute<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                             size_t nRow, size_t nCol, size_t nSec,
                             ePermute permute);
template
void array_permute<float   >(const float*    const arraySrc, float*    const arrayDst,
                             size_t nRow, size_t nCol, size_t nSec,
                             ePermute permute);
template
void array_permute<double  >(const double*   const arraySrc, double*   const arrayDst,
                             size_t nRow, size_t nCol, size_t nSec,
                             ePermute permute);
template
void array_permute<std::complex<float > >(const std::complex<float >* const arraySrc,
                                                std::complex<float >* const arrayDst,
                                                size_t nRow, size_t nCol, size_t nSec,
                                                ePermute permute);
template
void array_permute<std::complex<double> >(const std::complex<double>* const arraySrc,
                                                std::complex<double>* const arrayDst,
                                                size_t nRow, size_t nCol, size_t nSec,
                                                ePermute permute);

} // namespace gem

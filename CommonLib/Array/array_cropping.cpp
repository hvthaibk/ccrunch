/***********************************************************************
 *  File:       array_cropping.cpp
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
 * Cropping
 ****************************************/

// 1D
#ifndef NDEBUG
template <typename T1, typename T2>
void array_crop(const T1* const arraySrc, size_t nRowSrc,
                      T2* const arrayDst, size_t nRowDst,
                                          size_t nRowOff)
#else
template <typename T1, typename T2>
void array_crop(const T1* const arraySrc, size_t ,
                      T2* const arrayDst, size_t nRowDst,
                                          size_t nRowOff)
#endif
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0);
    assert(nRowDst > 0);
    assert(nRowDst + nRowOff <= nRowSrc);

    size_t    iRowDst;

    #pragma omp parallel for private(iRowDst)
    for (iRowDst = 0; iRowDst < nRowDst; iRowDst++) {
        arrayDst[iRowDst] = (T2) arraySrc[nRowOff+iRowDst];
    }
}

// instantiation
template
void array_crop<int32_t ,int32_t >(const int32_t*  const arraySrc, size_t nRowSrc,
                                         int32_t*  const arrayDst, size_t nRowDst,
                                                                   size_t nRowOff);
template
void array_crop<uint32_t,uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc,
                                         uint32_t* const arrayDst, size_t nRowDst,
                                                                   size_t nRowOff);
template
void array_crop<int64_t ,int64_t >(const int64_t*  const arraySrc, size_t nRowSrc,
                                         int64_t*  const arrayDst, size_t nRowDst,
                                                                   size_t nRowOff);
template
void array_crop<uint64_t,uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc,
                                         uint64_t* const arrayDst, size_t nRowDst,
                                                                   size_t nRowOff);
template
void array_crop<float   ,float   >(const float*    const arraySrc, size_t nRowSrc,
                                         float*    const arrayDst, size_t nRowDst,
                                                                   size_t nRowOff);
template
void array_crop<double  ,double  >(const double*   const arraySrc, size_t nRowSrc,
                                         double*   const arrayDst, size_t nRowDst,
                                                                   size_t nRowOff);

// 2D
#ifndef NDEBUG
template <typename T1, typename T2>
void array_crop(const T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                      T2* const arrayDst, size_t nRowDst, size_t nColDst,
                                          size_t nRowOff, size_t nColOff)
#else
template <typename T1, typename T2>
void array_crop(const T1* const arraySrc, size_t ,        size_t nColSrc,
                      T2* const arrayDst, size_t nRowDst, size_t nColDst,
                                          size_t nRowOff, size_t nColOff)
#endif
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0 && nColSrc > 0);
    assert(nRowDst > 0 && nColDst > 0);
    assert(nRowDst + nRowOff <= nRowSrc);
    assert(nColDst + nColOff <= nColSrc);

    size_t    iSrc1, iRowDst, iColDst;
    size_t    iDst1;
    size_t    nOff = nRowOff*nColSrc + nColOff;

    #pragma omp parallel for private(iRowDst,iSrc1,iDst1,iColDst)
    for (iRowDst = 0; iRowDst < nRowDst; iRowDst++) {
        iSrc1 = iRowDst*nColSrc + nOff;
        iDst1 = iRowDst*nColDst;

        for (iColDst = 0; iColDst < nColDst; iColDst++) {
            arrayDst[iDst1+iColDst] = (T2) arraySrc[iSrc1+iColDst];
        }
    }
}

// instantiation
template
void array_crop<int32_t ,int32_t >(const int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                         int32_t*  const arrayDst, size_t nRowDst, size_t nColDst,
                                                                   size_t nRowOff, size_t nColOff);
template
void array_crop<uint32_t,uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                         uint32_t* const arrayDst, size_t nRowDst, size_t nColDst,
                                                                   size_t nRowOff, size_t nColOff);
template
void array_crop<int64_t ,int64_t >(const int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                         int64_t*  const arrayDst, size_t nRowDst, size_t nColDst,
                                                                   size_t nRowOff, size_t nColOff);
template
void array_crop<uint64_t,uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                         uint64_t* const arrayDst, size_t nRowDst, size_t nColDst,
                                                                   size_t nRowOff, size_t nColOff);
template
void array_crop<float   ,float   >(const float*    const arraySrc, size_t nRowSrc, size_t nColSrc,
                                         float*    const arrayDst, size_t nRowDst, size_t nColDst,
                                                                   size_t nRowOff, size_t nColOff);
template
void array_crop<double  ,double  >(const double*   const arraySrc, size_t nRowSrc, size_t nColSrc,
                                         double*   const arrayDst, size_t nRowDst, size_t nColDst,
                                                                   size_t nRowOff, size_t nColOff);

// 3D
#ifndef NDEBUG
template <typename T1, typename T2>
void array_crop(const T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                      T2* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                          size_t nRowOff, size_t nColOff, size_t nSecOff)
#else
template <typename T1, typename T2>
void array_crop(const T1* const arraySrc, size_t ,        size_t nColSrc, size_t nSecSrc,
                      T2* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                          size_t nRowOff, size_t nColOff, size_t nSecOff)
#endif
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);
    assert(nRowDst > 0 && nColDst > 0 && nSecDst > 0);
    assert(nRowDst + nRowOff <= nRowSrc);
    assert(nColDst + nColOff <= nColSrc);
    assert(nSecDst + nSecOff <= nSecSrc);

    size_t    iSrc1, iSrc2, iRowDst, iColDst, iSecDst;
    size_t    iDst1, iDst2;
    size_t    nColSecSrc = nColSrc*nSecSrc;
    size_t    nColSecDst = nColDst*nSecDst;
    size_t    nOff = (nRowOff*nColSrc + nColOff)*nSecSrc + nSecOff;

    #pragma omp parallel for private(iRowDst,iSrc1,iDst1,iColDst,iSrc2,iDst2,iSecDst)
    for (iRowDst = 0; iRowDst < nRowDst; iRowDst++) {
        iSrc1 = iRowDst*nColSecSrc + nOff;
        iDst1 = iRowDst*nColSecDst;

        for (iColDst = 0; iColDst < nColDst; iColDst++) {
            iSrc2 = iSrc1 + iColDst*nSecSrc;
            iDst2 = iDst1 + iColDst*nSecDst;

            for (iSecDst = 0; iSecDst < nSecDst; iSecDst++) {
                arrayDst[iDst2+iSecDst] = (T2) arraySrc[iSrc2+iSecDst];
            }
        }
    }
}

// instantiation
template
void array_crop<int32_t ,int32_t >(const int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                         int32_t*  const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                   size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_crop<uint32_t,uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                         uint32_t* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                   size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_crop<int64_t ,int64_t >(const int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                         int64_t*  const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                   size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_crop<uint64_t,uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                         uint64_t* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                   size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_crop<float   ,float   >(const float*    const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                         float*    const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                   size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_crop<double  ,double  >(const double*   const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                         double*   const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                   size_t nRowOff, size_t nColOff, size_t nSecOff);

} // namespace gem

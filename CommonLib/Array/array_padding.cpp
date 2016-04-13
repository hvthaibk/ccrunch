/***********************************************************************
 *  File:       array_padding.cpp
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
 * Padding
 ****************************************/

// 1D
template <typename T1, typename T2>
void array_pad(const T1* const arraySrc, size_t nRowSrc,
                     T2* const arrayDst, size_t nRowDst,
                                         size_t nRowOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0);
    assert(nRowDst > 0);
    assert(nRowSrc + nRowOff <= nRowDst);

    size_t    iRowSrc;

    array_memset(arrayDst, 0, nRowDst);

    #pragma omp parallel for private(iRowSrc)
    for (iRowSrc = 0; iRowSrc < nRowSrc; iRowSrc++) {
        arrayDst[nRowOff+iRowSrc] = (T2) arraySrc[iRowSrc];
    }
}

// instantiation
template
void array_pad<int32_t ,int32_t >(const int32_t*  const arraySrc, size_t nRowSrc,
                                        int32_t*  const arrayDst, size_t nRowDst,
                                                                  size_t nRowOff);
template
void array_pad<uint32_t,uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc,
                                        uint32_t* const arrayDst, size_t nRowDst,
                                                                  size_t nRowOff);
template
void array_pad<int64_t ,int64_t >(const int64_t*  const arraySrc, size_t nRowSrc,
                                        int64_t*  const arrayDst, size_t nRowDst,
                                                                  size_t nRowOff);
template
void array_pad<uint64_t,uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc,
                                        uint64_t* const arrayDst, size_t nRowDst,
                                                                  size_t nRowOff);
template
void array_pad<float   ,float   >(const float*    const arraySrc, size_t nRowSrc,
                                        float*    const arrayDst, size_t nRowDst,
                                                                  size_t nRowOff);
template
void array_pad<double  ,double  >(const double*   const arraySrc, size_t nRowSrc,
                                        double*   const arrayDst, size_t nRowDst,
                                                                  size_t nRowOff);
template
void array_pad<std::complex<float >,std::complex<float > >
        (const std::complex<float >* const arraySrc, size_t nRowSrc,
               std::complex<float >* const arrayDst, size_t nRowDst,
                                                     size_t nRowOff);
template
void array_pad<std::complex<double>,std::complex<double> >
        (const std::complex<double>* const arraySrc, size_t nRowSrc,
               std::complex<double>* const arrayDst, size_t nRowDst,
                                                     size_t nRowOff);

// 2D
template <typename T1, typename T2>
void array_pad(const T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                     T2* const arrayDst, size_t nRowDst, size_t nColDst,
                                         size_t nRowOff, size_t nColOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0 && nColSrc > 0);
    assert(nRowDst > 0 && nColDst > 0);
    assert(nRowSrc + nRowOff <= nRowDst);
    assert(nColSrc + nColOff <= nColDst);

    size_t    iSrc1, iRowSrc, iColSrc;
    size_t    iDst1;
    size_t    nOff = nRowOff*nColDst + nColOff;

    array_memset(arrayDst, 0, nRowDst*nColDst);

    #pragma omp parallel for private(iRowSrc,iSrc1,iDst1,iColSrc)
    for (iRowSrc = 0; iRowSrc < nRowSrc; iRowSrc++) {
        iSrc1 = iRowSrc*nColSrc;
        iDst1 = iRowSrc*nColDst + nOff;

        for (iColSrc = 0; iColSrc < nColSrc; iColSrc++) {
            arrayDst[iDst1+iColSrc] = (T2) arraySrc[iSrc1+iColSrc];
        }
    }
}

// instantiation
template
void array_pad<int32_t ,int32_t >(const int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                        int32_t*  const arrayDst, size_t nRowDst, size_t nColDst,
                                                                  size_t nRowOff, size_t nColOff);
template
void array_pad<uint32_t,uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                        uint32_t* const arrayDst, size_t nRowDst, size_t nColDst,
                                                                  size_t nRowOff, size_t nColOff);
template
void array_pad<int64_t ,int64_t >(const int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                        int64_t*  const arrayDst, size_t nRowDst, size_t nColDst,
                                                                  size_t nRowOff, size_t nColOff);
template
void array_pad<uint64_t,uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                        uint64_t* const arrayDst, size_t nRowDst, size_t nColDst,
                                                                  size_t nRowOff, size_t nColOff);
template
void array_pad<float   ,float   >(const float*    const arraySrc, size_t nRowSrc, size_t nColSrc,
                                        float*    const arrayDst, size_t nRowDst, size_t nColDst,
                                                                  size_t nRowOff, size_t nColOff);
template
void array_pad<double  ,double  >(const double*   const arraySrc, size_t nRowSrc, size_t nColSrc,
                                        double*   const arrayDst, size_t nRowDst, size_t nColDst,
                                                                  size_t nRowOff, size_t nColOff);
template
void array_pad<std::complex<float >,std::complex<float > >
        (const std::complex<float >* const arraySrc, size_t nRowSrc, size_t nColSrc,
               std::complex<float >* const arrayDst, size_t nRowDst, size_t nColDst,
                                                     size_t nRowOff, size_t nColOff);
template
void array_pad<std::complex<double>,std::complex<double> >
        (const std::complex<double>* const arraySrc, size_t nRowSrc, size_t nColSrc,
               std::complex<double>* const arrayDst, size_t nRowDst, size_t nColDst,
                                                     size_t nRowOff, size_t nColOff);

// 3D
template <typename T1, typename T2>
void array_pad(const T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                     T2* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                         size_t nRowOff, size_t nColOff, size_t nSecOff)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);
    assert(nRowDst > 0 && nColDst > 0 && nSecDst > 0);
    assert(nRowSrc + nRowOff <= nRowDst);
    assert(nColSrc + nColOff <= nColDst);
    assert(nSecSrc + nSecOff <= nSecDst);

    size_t    iSrc1, iSrc2, iRowSrc, iColSrc, iSecSrc;
    size_t    iDst1, iDst2;
    size_t    nColSecSrc = nColSrc*nSecSrc;
    size_t    nColSecDst = nColDst*nSecDst;
    size_t    nOff = (nRowOff*nColDst + nColOff)*nSecDst + nSecOff;

    array_memset(arrayDst, 0, nRowDst*nColDst*nSecDst);

    #pragma omp parallel for private(iRowSrc,iSrc1,iDst1,iColSrc,iSrc2,iDst2,iSecSrc)
    for (iRowSrc = 0; iRowSrc < nRowSrc; iRowSrc++) {
        iSrc1 = iRowSrc*nColSecSrc;
        iDst1 = iRowSrc*nColSecDst + nOff;

        for (iColSrc = 0; iColSrc < nColSrc; iColSrc++) {
            iSrc2 = iSrc1 + iColSrc*nSecSrc;
            iDst2 = iDst1 + iColSrc*nSecDst;

            for (iSecSrc = 0; iSecSrc < nSecSrc; iSecSrc++) {
                arrayDst[iDst2+iSecSrc] = (T2) arraySrc[iSrc2+iSecSrc];
            }
        }
    }
}

// instantiation
template
void array_pad<int32_t ,int32_t >(const int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                        int32_t*  const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                  size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_pad<uint32_t,uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                        uint32_t* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                  size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_pad<int64_t ,int64_t >(const int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                        int64_t*  const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                  size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_pad<uint64_t,uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                        uint64_t* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                  size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_pad<float   ,float   >(const float*    const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                        float*    const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                  size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_pad<double  ,double  >(const double*   const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                        double*   const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                                  size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_pad<std::complex<float >,std::complex<float > >
        (const std::complex<float >* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
               std::complex<float >* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                     size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_pad<std::complex<double>,std::complex<double> >
        (const std::complex<double>* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
               std::complex<double>* const arrayDst, size_t nRowDst, size_t nColDst, size_t nSecDst,
                                                     size_t nRowOff, size_t nColOff, size_t nSecOff);

} // namespace gem

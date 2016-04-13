/***********************************************************************
 *  File:       array_replacement.cpp
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
 * Replacement
 ****************************************/

// 1D
#ifndef NDEBUG
template <typename T1, typename T2>
void array_replace(      T1* const arraySrc, size_t nRowSrc,
                   const T2* const arrayRep, size_t nRowRep,
                                             size_t nRowOff)
#else
template <typename T1, typename T2>
void array_replace(      T1* const arraySrc, size_t ,
                   const T2* const arrayRep, size_t nRowRep,
                                             size_t nRowOff)
#endif
{
    assert(arraySrc != NULL);
    assert(arrayRep != NULL);
    assert(nRowSrc > 0);
    assert(nRowRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);

    size_t    iRowRep;

    #pragma omp parallel for
    for (iRowRep = 0; iRowRep < nRowRep; iRowRep++) {
        arraySrc[iRowRep+nRowOff] = (T1) arrayRep[iRowRep];
    }
}

// instantiation
template
void array_replace<int32_t ,int32_t >(      int32_t*  const arraySrc, size_t nRowSrc,
                                      const int32_t*  const arrayRep, size_t nRowRep,
                                                                      size_t nRowOff);
template
void array_replace<uint32_t,uint32_t>(      uint32_t* const arraySrc, size_t nRowSrc,
                                      const uint32_t* const arrayRep, size_t nRowRep,
                                                                      size_t nRowOff);
template
void array_replace<int64_t ,int64_t >(      int64_t*  const arraySrc, size_t nRowSrc,
                                      const int64_t*  const arrayRep, size_t nRowRep,
                                                                      size_t nRowOff);
template
void array_replace<uint64_t,uint64_t>(      uint64_t* const arraySrc, size_t nRowSrc,
                                      const uint64_t* const arrayRep, size_t nRowRep,
                                                                      size_t nRowOff);
template
void array_replace<float   ,float   >(      float*    const arraySrc, size_t nRowSrc,
                                      const float*    const arrayRep, size_t nRowRep,
                                                                      size_t nRowOff);
template
void array_replace<double  ,double  >(      double*   const arraySrc, size_t nRowSrc,
                                      const double*   const arrayRep, size_t nRowRep,
                                                                      size_t nRowOff);

// 1D
#ifndef NDEBUG
template <typename T1, typename T2>
void array_replace(T1* const arraySrc, size_t nRowSrc,
                   T2        value,    size_t nRowRep,
                                       size_t nRowOff)
#else
template <typename T1, typename T2>
void array_replace(T1* const arraySrc, size_t ,
                   T2        value,    size_t nRowRep,
                                       size_t nRowOff)
#endif
{
    assert(arraySrc != NULL);
    assert(nRowSrc > 0);
    assert(nRowRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);

    size_t    iRowRep;

    #pragma omp parallel for
    for (iRowRep = 0; iRowRep < nRowRep; iRowRep++) {
        arraySrc[iRowRep+nRowOff] = (T1) value;
    }
}

// instantiation
template
void array_replace<int32_t ,int32_t >(int32_t*  const arraySrc, size_t nRowSrc,
                                      int32_t         value,    size_t nRowRep,
                                                                size_t nRowOff);
template
void array_replace<uint32_t,uint32_t>(uint32_t* const arraySrc, size_t nRowSrc,
                                      uint32_t        value,    size_t nRowRep,
                                                                size_t nRowOff);
template
void array_replace<int64_t ,int64_t >(int64_t*  const arraySrc, size_t nRowSrc,
                                      int64_t         value,    size_t nRowRep,
                                                                size_t nRowOff);
template
void array_replace<uint64_t,uint64_t>(uint64_t* const arraySrc, size_t nRowSrc,
                                      uint64_t        value,    size_t nRowRep,
                                                                size_t nRowOff);
template
void array_replace<float   ,float   >(float*    const arraySrc, size_t nRowSrc,
                                      float           value,    size_t nRowRep,
                                                                size_t nRowOff);
template
void array_replace<double  ,double  >(double*   const arraySrc, size_t nRowSrc,
                                      double          value,    size_t nRowRep,
                                                                size_t nRowOff);

// 2D
#ifndef NDEBUG
template <typename T1, typename T2>
void array_replace(      T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                   const T2* const arrayRep, size_t nRowRep, size_t nColRep,
                                             size_t nRowOff, size_t nColOff)
#else
template <typename T1, typename T2>
void array_replace(      T1* const arraySrc, size_t        , size_t nColSrc,
                   const T2* const arrayRep, size_t nRowRep, size_t nColRep,
                                             size_t nRowOff, size_t nColOff)
#endif
{
    assert(arraySrc != NULL);
    assert(arrayRep != NULL);
    assert(nRowSrc > 0 && nColSrc > 0);
    assert(nRowRep > 0 && nColRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);
    assert(nColRep + nColOff <= nColSrc);

    size_t    iRowRep, iColRep;
    size_t    iSrc1;
    size_t    iRep1;
    size_t    nOff = nRowOff*nColSrc + nColOff;

    #pragma omp parallel for private(iRowRep,iSrc1,iRep1,iColRep)
    for (iRowRep = 0; iRowRep < nRowRep; iRowRep++) {
        iSrc1 = iRowRep*nColSrc + nOff;
        iRep1 = iRowRep*nColRep;

        for (iColRep = 0; iColRep < nColRep; iColRep++) {
            arraySrc[iSrc1+iColRep] = (T1) arrayRep[iRep1+iColRep];
        }
    }
}

// instantiation
template
void array_replace<int32_t ,int32_t >(      int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      const int32_t*  const arrayRep, size_t nRowRep, size_t nColRep,
                                                                      size_t nRowOff, size_t nColOff);
template
void array_replace<uint32_t,uint32_t>(      uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      const uint32_t* const arrayRep, size_t nRowRep, size_t nColRep,
                                                                      size_t nRowOff, size_t nColOff);
template
void array_replace<int64_t ,int64_t >(      int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      const int64_t*  const arrayRep, size_t nRowRep, size_t nColRep,
                                                                      size_t nRowOff, size_t nColOff);
template
void array_replace<uint64_t,uint64_t>(      uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      const uint64_t* const arrayRep, size_t nRowRep, size_t nColRep,
                                                                      size_t nRowOff, size_t nColOff);
template
void array_replace<float   ,float   >(      float*    const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      const float*    const arrayRep, size_t nRowRep, size_t nColRep,
                                                                      size_t nRowOff, size_t nColOff);
template
void array_replace<double  ,double  >(      double*   const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      const double*   const arrayRep, size_t nRowRep, size_t nColRep,
                                                                      size_t nRowOff, size_t nColOff);

// 2D
#ifndef NDEBUG
template <typename T1, typename T2>
void array_replace(T1* const arraySrc, size_t nRowSrc, size_t nColSrc,
                   T2        value,    size_t nRowRep, size_t nColRep,
                                       size_t nRowOff, size_t nColOff)
#else
template <typename T1, typename T2>
void array_replace(T1* const arraySrc, size_t        , size_t nColSrc,
                   T2        value,    size_t nRowRep, size_t nColRep,
                                       size_t nRowOff, size_t nColOff)
#endif
{
    assert(arraySrc != NULL);
    assert(nRowSrc > 0 && nColSrc > 0);
    assert(nRowRep > 0 && nColRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);
    assert(nColRep + nColOff <= nColSrc);

    size_t    iRowRep, iColRep;
    size_t    iSrc1;
    size_t    nOff = nRowOff*nColSrc + nColOff;

    #pragma omp parallel for private(iRowRep,iSrc1,iColRep)
    for (iRowRep = 0; iRowRep < nRowRep; iRowRep++) {
        iSrc1 = iRowRep*nColSrc + nOff;

        for (iColRep = 0; iColRep < nColRep; iColRep++) {
            arraySrc[iSrc1+iColRep] = (T1) value;
        }
    }
}

// instantiation
template
void array_replace<int32_t ,int32_t >(int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      int32_t         value,    size_t nRowRep, size_t nColRep,
                                                                size_t nRowOff, size_t nColOff);
template
void array_replace<uint32_t,uint32_t>(uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      uint32_t        value,    size_t nRowRep, size_t nColRep,
                                                                size_t nRowOff, size_t nColOff);
template
void array_replace<int64_t ,int64_t >(int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      int64_t         value,    size_t nRowRep, size_t nColRep,
                                                                size_t nRowOff, size_t nColOff);
template
void array_replace<uint64_t,uint64_t>(uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      uint64_t        value,    size_t nRowRep, size_t nColRep,
                                                                size_t nRowOff, size_t nColOff);
template
void array_replace<float   ,float   >(float*    const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      float           value,    size_t nRowRep, size_t nColRep,
                                                                size_t nRowOff, size_t nColOff);
template
void array_replace<double  ,double  >(double*   const arraySrc, size_t nRowSrc, size_t nColSrc,
                                      double          value,    size_t nRowRep, size_t nColRep,
                                                                size_t nRowOff, size_t nColOff);

// 3D
#ifndef NDEBUG
template <typename T1, typename T2>
void array_replace(      T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                   const T2* const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                             size_t nRowOff, size_t nColOff, size_t nSecOff)
#else
template <typename T1, typename T2>
void array_replace(      T1* const arraySrc, size_t        , size_t nColSrc, size_t nSecSrc,
                   const T2* const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                             size_t nRowOff, size_t nColOff, size_t nSecOff)
#endif
{
    assert(arraySrc != NULL);
    assert(arrayRep != NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);
    assert(nRowRep > 0 && nColRep > 0 && nSecRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);
    assert(nColRep + nColOff <= nColSrc);
    assert(nSecRep + nSecOff <= nSecSrc);

    size_t    iRowRep, iColRep, iSecRep;
    size_t    iSrc1, iSrc2;
    size_t    iRep1, iRep2;
    size_t    nColSecSrc = nColSrc*nSecSrc;
    size_t    nColSecRep = nColRep*nSecRep;
    size_t    nOff = (nRowOff*nColSrc + nColOff)*nSecSrc + nSecOff;

    #pragma omp parallel for private(iRowRep,iSrc1,iRep1,iColRep,iSrc2,iRep2,iSecRep)
    for (iRowRep = 0; iRowRep < nRowRep; iRowRep++) {
        iSrc1 = iRowRep*nColSecSrc + nOff;
        iRep1 = iRowRep*nColSecRep;

        for (iColRep = 0; iColRep < nColRep; iColRep++) {
            iSrc2 = iSrc1 + iColRep*nSecSrc;
            iRep2 = iRep1 + iColRep*nSecRep;

            for (iSecRep = 0; iSecRep < nSecRep; iSecRep++) {
                arraySrc[iSrc2+iSecRep] = (T1) arrayRep[iRep2+iSecRep];
            }
        }
    }
}

// instantiation
template
void array_replace<int32_t ,int32_t >(      int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      const int32_t*  const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                      size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_replace<uint32_t,uint32_t>(      uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      const uint32_t* const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                      size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_replace<int64_t ,int64_t >(      int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      const int64_t*  const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                      size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_replace<uint64_t,uint64_t>(      uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      const uint64_t* const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                      size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_replace<float   ,float   >(      float*    const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      const float*    const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                      size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_replace<double  ,double  >(      double*   const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      const double*   const arrayRep, size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                      size_t nRowOff, size_t nColOff, size_t nSecOff);

// 3D
#ifndef NDEBUG
template <typename T1, typename T2>
void array_replace(T1* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                   T2        value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                       size_t nRowOff, size_t nColOff, size_t nSecOff)
#else
template <typename T1, typename T2>
void array_replace(T1* const arraySrc, size_t        , size_t nColSrc, size_t nSecSrc,
                   T2        value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                       size_t nRowOff, size_t nColOff, size_t nSecOff)
#endif
{
    assert(arraySrc != NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);
    assert(nRowRep > 0 && nColRep > 0 && nSecRep > 0);
    assert(nRowRep + nRowOff <= nRowSrc);
    assert(nColRep + nColOff <= nColSrc);
    assert(nSecRep + nSecOff <= nSecSrc);

    size_t    iRowRep, iColRep, iSecRep;
    size_t    iSrc1, iSrc2;
    size_t    nColSecSrc = nColSrc*nSecSrc;
    size_t    nOff = (nRowOff*nColSrc + nColOff)*nSecSrc + nSecOff;

    #pragma omp parallel for private(iRowRep,iSrc1,iColRep,iSrc2,iSecRep)
    for (iRowRep = 0; iRowRep < nRowRep; iRowRep++) {
        iSrc1 = iRowRep*nColSecSrc + nOff;

        for (iColRep = 0; iColRep < nColRep; iColRep++) {
            iSrc2 = iSrc1 + iColRep*nSecSrc;

            for (iSecRep = 0; iSecRep < nSecRep; iSecRep++) {
                arraySrc[iSrc2+iSecRep] = (T1) value;
            }
        }
    }
}

// instantiation
template
void array_replace<int32_t ,int32_t >(int32_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      int32_t         value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_replace<uint32_t,uint32_t>(uint32_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      uint32_t        value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_replace<int64_t ,int64_t >(int64_t*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      int64_t         value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_replace<uint64_t,uint64_t>(uint64_t* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      uint64_t        value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_replace<float   ,float   >(float*    const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      float           value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                size_t nRowOff, size_t nColOff, size_t nSecOff);
template
void array_replace<double  ,double  >(double*   const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                      double          value,    size_t nRowRep, size_t nColRep, size_t nSecRep,
                                                                size_t nRowOff, size_t nColOff, size_t nSecOff);

} // namespace gem

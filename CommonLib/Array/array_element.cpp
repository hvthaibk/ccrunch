/***********************************************************************
 *  File:       array_element.cpp
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
 * Element-wise operations
 ****************************************/

// find
template <typename T>
void array_find(const T*      const arraySrc, size_t nRow,
                      size_t* const arrayDst,
                      size_t&       nFound)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    size_t    iRowSrc, iRowDst = 0;

    for (iRowSrc = 0; iRowSrc < nRow; iRowSrc++) {
        if (arraySrc[iRowSrc] > 0) {
            arrayDst[iRowDst++] = iRowSrc;
        }
    }

    nFound = iRowDst;
}

// instantiation
template
void array_find<int32_t >(const int32_t*  const arraySrc, size_t nRow,
                                size_t*   const arrayDst,
                                size_t&         nFound);
template
void array_find<uint32_t>(const uint32_t* const arraySrc, size_t nRow,
                                size_t*   const arrayDst,
                                size_t&         nFound);
template
void array_find<int64_t >(const int64_t*  const arraySrc, size_t nRow,
                                size_t*   const arrayDst,
                                size_t&         nFound);
template
void array_find<uint64_t>(const uint64_t* const arraySrc, size_t nRow,
                                size_t*   const arrayDst,
                                size_t&         nFound);
template
void array_find<float   >(const float*    const arraySrc, size_t nRow,
                                size_t*   const arrayDst,
                                size_t&         nFound);
template
void array_find<double  >(const double*   const arraySrc, size_t nRow,
                                size_t*   const arrayDst,
                                size_t&         nFound);

// ismember
template <typename T>
void array_ismember(const T* const arraySrc, size_t nRowSrc,
                    const T* const arraySet, size_t nRowSet,
                          T* const arrayMem)
{
    assert(arraySrc != NULL);
    assert(arraySet != NULL);
    assert(arrayMem != NULL);
    assert(nRowSrc > 0);
    assert(nRowSet > 0);

    #pragma omp parallel for
    for (size_t iRowSrc = 0; iRowSrc < nRowSrc; iRowSrc++) {
        arrayMem[iRowSrc] = (T) 0;

        for (size_t iRowSet = 0; iRowSet < nRowSet; iRowSet++) {
            if (arraySrc[iRowSrc] == arraySet[iRowSet]) {
                arrayMem[iRowSrc] = (T) 1;
                break;
            }
        }
    }
}

// instantiation
template
void array_ismember<int32_t >(const int32_t*  const arraySrc, size_t nRowSrc,
                              const int32_t*  const arraySet, size_t nRowSet,
                                    int32_t*  const arrayMem);
template
void array_ismember<uint32_t>(const uint32_t* const arraySrc, size_t nRowSrc,
                              const uint32_t* const arraySet, size_t nRowSet,
                                    uint32_t* const arrayMem);
template
void array_ismember<int64_t >(const int64_t*  const arraySrc, size_t nRowSrc,
                              const int64_t*  const arraySet, size_t nRowSet,
                                    int64_t*  const arrayMem);
template
void array_ismember<uint64_t>(const uint64_t* const arraySrc, size_t nRowSrc,
                              const uint64_t* const arraySet, size_t nRowSet,
                                    uint64_t* const arrayMem);
template
void array_ismember<float   >(const float*    const arraySrc, size_t nRowSrc,
                              const float*    const arraySet, size_t nRowSet,
                                    float*    const arrayMem);
template
void array_ismember<double  >(const double*   const arraySrc, size_t nRowSrc,
                              const double*   const arraySet, size_t nRowSet,
                                    double*   const arrayMem);

// set value
template <typename T>
void array_setval(T* const array, T value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] = value;
    }
}

// instantiation
template
void array_setval<int32_t >(int32_t*  const array, int32_t  value, size_t nRow);
template
void array_setval<uint32_t>(uint32_t* const array, uint32_t value, size_t nRow);
template
void array_setval<int64_t >(int64_t*  const array, int64_t  value, size_t nRow);
template
void array_setval<uint64_t>(uint64_t* const array, uint64_t value, size_t nRow);
template
void array_setval<float   >(float*    const array, float    value, size_t nRow);
template
void array_setval<double  >(double*   const array, double   value, size_t nRow);

// type cast
template <typename T1, typename T2>
void array_typecast(const T1* const arraySrc, size_t nRow,
                          T2* const arrayDst)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(arraySrc != (reinterpret_cast<T1*>(arrayDst)));
    assert(nRow > 0);

    size_t    iRow;

    #pragma omp parallel for
    for (iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = (T2) arraySrc[iRow];
    }
}

// instantiation
template
void array_typecast<uint8_t ,int8_t  >(const uint8_t*  const arraySrc, size_t nRow, int8_t*   const arrayDst);
template
void array_typecast<uint8_t ,uint8_t >(const uint8_t*  const arraySrc, size_t nRow, uint8_t*  const arrayDst);
template
void array_typecast<uint8_t ,int16_t >(const uint8_t*  const arraySrc, size_t nRow, int16_t*  const arrayDst);
template
void array_typecast<uint8_t ,uint16_t>(const uint8_t*  const arraySrc, size_t nRow, uint16_t* const arrayDst);
template
void array_typecast<uint8_t ,int32_t >(const uint8_t*  const arraySrc, size_t nRow, int32_t*  const arrayDst);
template
void array_typecast<uint8_t ,uint32_t>(const uint8_t*  const arraySrc, size_t nRow, uint32_t* const arrayDst);
template
void array_typecast<uint8_t ,int64_t >(const uint8_t*  const arraySrc, size_t nRow, int64_t*  const arrayDst);
template
void array_typecast<uint8_t ,uint64_t>(const uint8_t*  const arraySrc, size_t nRow, uint64_t* const arrayDst);
template
void array_typecast<uint8_t ,float   >(const uint8_t*  const arraySrc, size_t nRow, float*    const arrayDst);
template
void array_typecast<uint8_t ,double  >(const uint8_t*  const arraySrc, size_t nRow, double*   const arrayDst);
template
void array_typecast<uint16_t,int8_t  >(const uint16_t* const arraySrc, size_t nRow, int8_t*   const arrayDst);
template
void array_typecast<uint16_t,uint8_t >(const uint16_t* const arraySrc, size_t nRow, uint8_t*  const arrayDst);
template
void array_typecast<uint16_t,int16_t >(const uint16_t* const arraySrc, size_t nRow, int16_t*  const arrayDst);
template
void array_typecast<uint16_t,uint16_t>(const uint16_t* const arraySrc, size_t nRow, uint16_t* const arrayDst);
template
void array_typecast<uint16_t,int32_t >(const uint16_t* const arraySrc, size_t nRow, int32_t*  const arrayDst);
template
void array_typecast<uint16_t,uint32_t>(const uint16_t* const arraySrc, size_t nRow, uint32_t* const arrayDst);
template
void array_typecast<uint16_t,int64_t >(const uint16_t* const arraySrc, size_t nRow, int64_t*  const arrayDst);
template
void array_typecast<uint16_t,uint64_t>(const uint16_t* const arraySrc, size_t nRow, uint64_t* const arrayDst);
template
void array_typecast<uint16_t,float   >(const uint16_t* const arraySrc, size_t nRow, float*    const arrayDst);
template
void array_typecast<uint16_t,double  >(const uint16_t* const arraySrc, size_t nRow, double*   const arrayDst);
template
void array_typecast<uint32_t,int8_t  >(const uint32_t* const arraySrc, size_t nRow, int8_t*   const arrayDst);
template
void array_typecast<uint32_t,uint8_t >(const uint32_t* const arraySrc, size_t nRow, uint8_t*  const arrayDst);
template
void array_typecast<uint32_t,int16_t >(const uint32_t* const arraySrc, size_t nRow, int16_t*  const arrayDst);
template
void array_typecast<uint32_t,uint16_t>(const uint32_t* const arraySrc, size_t nRow, uint16_t* const arrayDst);
template
void array_typecast<uint32_t,int32_t >(const uint32_t* const arraySrc, size_t nRow, int32_t*  const arrayDst);
template
void array_typecast<uint32_t,uint32_t>(const uint32_t* const arraySrc, size_t nRow, uint32_t* const arrayDst);
template
void array_typecast<uint32_t,int64_t >(const uint32_t* const arraySrc, size_t nRow, int64_t*  const arrayDst);
template
void array_typecast<uint32_t,uint64_t>(const uint32_t* const arraySrc, size_t nRow, uint64_t* const arrayDst);
template
void array_typecast<uint32_t,float   >(const uint32_t* const arraySrc, size_t nRow, float*    const arrayDst);
template
void array_typecast<uint32_t,double  >(const uint32_t* const arraySrc, size_t nRow, double*   const arrayDst);
template
void array_typecast<uint64_t,int8_t  >(const uint64_t* const arraySrc, size_t nRow, int8_t*   const arrayDst);
template
void array_typecast<uint64_t,uint8_t >(const uint64_t* const arraySrc, size_t nRow, uint8_t*  const arrayDst);
template
void array_typecast<uint64_t,int16_t >(const uint64_t* const arraySrc, size_t nRow, int16_t*  const arrayDst);
template
void array_typecast<uint64_t,uint16_t>(const uint64_t* const arraySrc, size_t nRow, uint16_t* const arrayDst);
template
void array_typecast<uint64_t,int32_t >(const uint64_t* const arraySrc, size_t nRow, int32_t*  const arrayDst);
template
void array_typecast<uint64_t,uint32_t>(const uint64_t* const arraySrc, size_t nRow, uint32_t* const arrayDst);
template
void array_typecast<uint64_t,int64_t >(const uint64_t* const arraySrc, size_t nRow, int64_t*  const arrayDst);
template
void array_typecast<uint64_t,uint64_t>(const uint64_t* const arraySrc, size_t nRow, uint64_t* const arrayDst);
template
void array_typecast<uint64_t,float   >(const uint64_t* const arraySrc, size_t nRow, float*    const arrayDst);
template
void array_typecast<uint64_t,double  >(const uint64_t* const arraySrc, size_t nRow, double*   const arrayDst);
template
void array_typecast<int8_t  ,int8_t  >(const int8_t*   const arraySrc, size_t nRow, int8_t*   const arrayDst);
template
void array_typecast<int8_t  ,uint8_t >(const int8_t*   const arraySrc, size_t nRow, uint8_t*  const arrayDst);
template
void array_typecast<int8_t  ,int16_t >(const int8_t*   const arraySrc, size_t nRow, int16_t*  const arrayDst);
template
void array_typecast<int8_t  ,uint16_t>(const int8_t*   const arraySrc, size_t nRow, uint16_t* const arrayDst);
template
void array_typecast<int8_t  ,int32_t >(const int8_t*   const arraySrc, size_t nRow, int32_t*  const arrayDst);
template
void array_typecast<int8_t  ,uint32_t>(const int8_t*   const arraySrc, size_t nRow, uint32_t* const arrayDst);
template
void array_typecast<int8_t  ,int64_t >(const int8_t*   const arraySrc, size_t nRow, int64_t*  const arrayDst);
template
void array_typecast<int8_t  ,uint64_t>(const int8_t*   const arraySrc, size_t nRow, uint64_t* const arrayDst);
template
void array_typecast<int8_t  ,float   >(const int8_t*   const arraySrc, size_t nRow, float*    const arrayDst);
template
void array_typecast<int8_t  ,double  >(const int8_t*   const arraySrc, size_t nRow, double*   const arrayDst);
template
void array_typecast<int16_t ,int8_t  >(const int16_t*  const arraySrc, size_t nRow, int8_t*   const arrayDst);
template
void array_typecast<int16_t ,uint8_t >(const int16_t*  const arraySrc, size_t nRow, uint8_t*  const arrayDst);
template
void array_typecast<int16_t ,int16_t >(const int16_t*  const arraySrc, size_t nRow, int16_t*  const arrayDst);
template
void array_typecast<int16_t ,uint16_t>(const int16_t*  const arraySrc, size_t nRow, uint16_t* const arrayDst);
template
void array_typecast<int16_t ,int32_t >(const int16_t*  const arraySrc, size_t nRow, int32_t*  const arrayDst);
template
void array_typecast<int16_t ,uint32_t>(const int16_t*  const arraySrc, size_t nRow, uint32_t* const arrayDst);
template
void array_typecast<int16_t ,int64_t >(const int16_t*  const arraySrc, size_t nRow, int64_t*  const arrayDst);
template
void array_typecast<int16_t ,uint64_t>(const int16_t*  const arraySrc, size_t nRow, uint64_t* const arrayDst);
template
void array_typecast<int16_t ,float   >(const int16_t*  const arraySrc, size_t nRow, float*    const arrayDst);
template
void array_typecast<int16_t ,double  >(const int16_t*  const arraySrc, size_t nRow, double*   const arrayDst);
template
void array_typecast<int32_t ,int8_t  >(const int32_t*  const arraySrc, size_t nRow, int8_t*   const arrayDst);
template
void array_typecast<int32_t ,uint8_t >(const int32_t*  const arraySrc, size_t nRow, uint8_t*  const arrayDst);
template
void array_typecast<int32_t ,int16_t >(const int32_t*  const arraySrc, size_t nRow, int16_t*  const arrayDst);
template
void array_typecast<int32_t ,uint16_t>(const int32_t*  const arraySrc, size_t nRow, uint16_t* const arrayDst);
template
void array_typecast<int32_t ,int32_t >(const int32_t*  const arraySrc, size_t nRow, int32_t*  const arrayDst);
template
void array_typecast<int32_t ,uint32_t>(const int32_t*  const arraySrc, size_t nRow, uint32_t* const arrayDst);
template
void array_typecast<int32_t ,int64_t >(const int32_t*  const arraySrc, size_t nRow, int64_t*  const arrayDst);
template
void array_typecast<int32_t ,uint64_t>(const int32_t*  const arraySrc, size_t nRow, uint64_t* const arrayDst);
template
void array_typecast<int32_t ,float   >(const int32_t*  const arraySrc, size_t nRow, float*    const arrayDst);
template
void array_typecast<int32_t ,double  >(const int32_t*  const arraySrc, size_t nRow, double*   const arrayDst);
template
void array_typecast<int64_t ,int8_t  >(const int64_t*  const arraySrc, size_t nRow, int8_t*   const arrayDst);
template
void array_typecast<int64_t ,uint8_t >(const int64_t*  const arraySrc, size_t nRow, uint8_t*  const arrayDst);
template
void array_typecast<int64_t ,int16_t >(const int64_t*  const arraySrc, size_t nRow, int16_t*  const arrayDst);
template
void array_typecast<int64_t ,uint16_t>(const int64_t*  const arraySrc, size_t nRow, uint16_t* const arrayDst);
template
void array_typecast<int64_t ,int32_t >(const int64_t*  const arraySrc, size_t nRow, int32_t*  const arrayDst);
template
void array_typecast<int64_t ,uint32_t>(const int64_t*  const arraySrc, size_t nRow, uint32_t* const arrayDst);
template
void array_typecast<int64_t ,int64_t >(const int64_t*  const arraySrc, size_t nRow, int64_t*  const arrayDst);
template
void array_typecast<int64_t ,uint64_t>(const int64_t*  const arraySrc, size_t nRow, uint64_t* const arrayDst);
template
void array_typecast<int64_t ,float   >(const int64_t*  const arraySrc, size_t nRow, float*    const arrayDst);
template
void array_typecast<int64_t ,double  >(const int64_t*  const arraySrc, size_t nRow, double*   const arrayDst);
template
void array_typecast<float   ,int8_t  >(const float*    const arraySrc, size_t nRow, int8_t*   const arrayDst);
template
void array_typecast<float   ,uint8_t >(const float*    const arraySrc, size_t nRow, uint8_t*  const arrayDst);
template
void array_typecast<float   ,int16_t >(const float*    const arraySrc, size_t nRow, int16_t*  const arrayDst);
template
void array_typecast<float   ,uint16_t>(const float*    const arraySrc, size_t nRow, uint16_t* const arrayDst);
template
void array_typecast<float   ,int32_t >(const float*    const arraySrc, size_t nRow, int32_t*  const arrayDst);
template
void array_typecast<float   ,uint32_t>(const float*    const arraySrc, size_t nRow, uint32_t* const arrayDst);
template
void array_typecast<float   ,int64_t >(const float*    const arraySrc, size_t nRow, int64_t*  const arrayDst);
template
void array_typecast<float   ,uint64_t>(const float*    const arraySrc, size_t nRow, uint64_t* const arrayDst);
template
void array_typecast<float   ,float   >(const float*    const arraySrc, size_t nRow, float*    const arrayDst);
template
void array_typecast<float   ,double  >(const float*    const arraySrc, size_t nRow, double*   const arrayDst);
template
void array_typecast<double  ,int8_t  >(const double*   const arraySrc, size_t nRow, int8_t*   const arrayDst);
template
void array_typecast<double  ,uint8_t >(const double*   const arraySrc, size_t nRow, uint8_t*  const arrayDst);
template
void array_typecast<double  ,int16_t >(const double*   const arraySrc, size_t nRow, int16_t*  const arrayDst);
template
void array_typecast<double  ,uint16_t>(const double*   const arraySrc, size_t nRow, uint16_t* const arrayDst);
template
void array_typecast<double  ,int32_t >(const double*   const arraySrc, size_t nRow, int32_t*  const arrayDst);
template
void array_typecast<double  ,uint32_t>(const double*   const arraySrc, size_t nRow, uint32_t* const arrayDst);
template
void array_typecast<double  ,int64_t >(const double*   const arraySrc, size_t nRow, int64_t*  const arrayDst);
template
void array_typecast<double  ,uint64_t>(const double*   const arraySrc, size_t nRow, uint64_t* const arrayDst);
template
void array_typecast<double  ,float   >(const double*   const arraySrc, size_t nRow, float*    const arrayDst);
template
void array_typecast<double  ,double  >(const double*   const arraySrc, size_t nRow, double*   const arrayDst);

// threshold elements to 0/1
template <typename T>
void array_threshold(T* const array, size_t nRow,
                     T        cutoff)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] = (T) ((array[iRow] >= cutoff) ? 1 : 0);
    }
}

template <typename T1, typename T2>
void array_threshold(const T1* const arraySrc, size_t nRow,
                           T1        cutoff,
                           T2* const arrayDst)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    size_t    iRow;

    #pragma omp parallel for
    for (iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = (T2) ((arraySrc[iRow] >= cutoff) ? 1 : 0);
    }
}

template <typename T1, typename T2>
void array_threshold(const T1* const arraySrc, size_t nRow,
                     const T1* const arrayThr,
                           T2* const arrayDst)
{
    assert(arraySrc != NULL);
    assert(arrayThr != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = (T2) ((arraySrc[iRow] >= arrayThr[iRow]) ? 1 : 0);
    }
}

// instantiation
template
void array_threshold<int32_t >(int32_t*  const array, size_t nRow,
                               int32_t         cutoff);
template
void array_threshold<uint32_t>(uint32_t* const array, size_t nRow,
                               uint32_t        cutoff);
template
void array_threshold<int64_t >(int64_t*  const array, size_t nRow,
                               int64_t         cutoff);
template
void array_threshold<uint64_t>(uint64_t* const array, size_t nRow,
                               uint64_t        cutoff);
template
void array_threshold<float   >(float*    const array, size_t nRow,
                               float           cutoff);
template
void array_threshold<double  >(double*   const array, size_t nRow,
                               double          cutoff);
template
void array_threshold<int32_t , int32_t >(const int32_t*  const arraySrc, size_t nRow,
                                               int32_t         cutoff,
                                               int32_t*  const arrayDst);
template
void array_threshold<uint32_t, uint32_t>(const uint32_t* const arraySrc, size_t nRow,
                                               uint32_t        cutoff,
                                               uint32_t* const arrayDst);
template
void array_threshold<int64_t , int64_t >(const int64_t*  const arraySrc, size_t nRow,
                                               int64_t         cutoff,
                                               int64_t*  const arrayDst);
template
void array_threshold<uint64_t, uint64_t>(const uint64_t* const arraySrc, size_t nRow,
                                               uint64_t        cutoff,
                                               uint64_t* const arrayDst);
template
void array_threshold<float   , float   >(const float*    const arraySrc, size_t nRow,
                                               float           cutoff,
                                               float*    const arrayDst);
template
void array_threshold<double  , double  >(const double*   const arraySrc, size_t nRow,
                                               double          cutoff,
                                               double*   const arrayDst);
template
void array_threshold<int32_t , int32_t >(const int32_t*  const arraySrc, size_t nRow,
                                         const int32_t*  const arrayThr,
                                               int32_t*  const arrayDst);
template
void array_threshold<uint32_t, uint32_t>(const uint32_t* const arraySrc, size_t nRow,
                                         const uint32_t* const arrayThr,
                                               uint32_t* const arrayDst);
template
void array_threshold<int64_t , int64_t >(const int64_t*  const arraySrc, size_t nRow,
                                         const int64_t*  const arrayThr,
                                               int64_t*  const arrayDst);
template
void array_threshold<uint64_t, uint64_t>(const uint64_t* const arraySrc, size_t nRow,
                                         const uint64_t* const arrayThr,
                                               uint64_t* const arrayDst);
template
void array_threshold<float   , float   >(const float*    const arraySrc, size_t nRow,
                                         const float*    const arrayThr,
                                               float*    const arrayDst);
template
void array_threshold<double  , double  >(const double*   const arraySrc, size_t nRow,
                                         const double*   const arrayThr,
                                               double*   const arrayDst);

// cutoff small elements to value
template <typename T>
void array_cutoff(T* const array, size_t nRow,
                  T        cutoff,
                  T        value)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (std::abs(array[iRow]) <= cutoff) {
            array[iRow] = value;
        }
    }
}

template <typename T>
void array_cutoff(const T* const arraySrc, size_t nRow,
                        T        cutoff,
                        T* const arrayDst,
                        T        value)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = std::abs(arraySrc[iRow]) <= cutoff ? value : arraySrc[iRow];
    }
}

template <typename T>
void array_cutoff(const T* const arraySrc, size_t nRow,
                  const T* const arrayThr,
                        T* const arrayDst,
                        T        value)
{
    assert(arraySrc != NULL);
    assert(arrayThr != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = std::abs(arraySrc[iRow]) <= arrayThr[iRow] ? value : arraySrc[iRow];
    }
}

// instantiation
template
void array_cutoff<float >(float*  const array, size_t nRow,
                          float         cutoff,
                          float         value);
template
void array_cutoff<double>(double* const array, size_t nRow,
                          double        cutoff,
                          double        value);
template
void array_cutoff<float >(const float*  const arraySrc, size_t nRow,
                                float         cutoff,
                                float*  const arrayDst,
                                float         value);
template
void array_cutoff<double>(const double* const arraySrc, size_t nRow,
                                double        cutoff,
                                double* const arrayDst,
                                double        value);
template
void array_cutoff<float >(const float*  const arraySrc, size_t nRow,
                          const float*  const arrayThr,
                                float*  const arrayDst,
                                float         value);
template
void array_cutoff<double>(const double* const arraySrc, size_t nRow,
                          const double* const arrayThr,
                                double* const arrayDst,
                                double        value);

// suppress big elements to value
template <typename T>
void array_suppress(T* const array, size_t nRow,
                    T        cutoff,
                    T        value)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (std::abs(array[iRow]) > cutoff) {
            array[iRow] = value;
        }
    }
}

template <typename T>
void array_suppress(const T* const arraySrc, size_t nRow,
                          T        cutoff,
                          T* const arrayDst,
                          T        value)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = std::abs(arraySrc[iRow]) > cutoff ? value : arraySrc[iRow];
    }
}

template <typename T>
void array_suppress(const T* const arraySrc, size_t nRow,
                    const T* const arrayThr,
                          T* const arrayDst,
                          T        value)
{
    assert(arraySrc != NULL);
    assert(arrayThr != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = std::abs(arraySrc[iRow]) > arrayThr[iRow] ? value : arraySrc[iRow];
    }
}

// instantiation
template
void array_suppress<float >(float*  const array, size_t nRow,
                            float         cutoff,
                            float         value);
template
void array_suppress<double>(double* const array, size_t nRow,
                            double        cutoff,
                            double        value);
template
void array_suppress<float >(const float*  const arraySrc, size_t nRow,
                                  float         cutoff,
                                  float*  const arrayDst,
                                  float         value);
template
void array_suppress<double>(const double* const arraySrc, size_t nRow,
                                  double        cutoff,
                                  double* const arrayDst,
                                  double        value);
template
void array_suppress<float >(const float*  const arraySrc, size_t nRow,
                            const float*  const arrayThr,
                                  float*  const arrayDst,
                                  float         value);
template
void array_suppress<double>(const double* const arraySrc, size_t nRow,
                            const double* const arrayThr,
                                  double* const arrayDst,
                                  double        value);

// remove zero elements by adding noise
template <typename T>
void array_remove_rezo(const T* const arrayNoise, size_t nRow,
                             T* const arrayData,
                             T        cutoff)
{
    assert(arrayNoise != NULL);
    assert(arrayData  != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (std::abs(arrayData[iRow]) <= cutoff) {
            arrayData[iRow] += arrayNoise[iRow];
        }
    }
}

// instantiation
template
void array_remove_rezo<float >(const float*  const arrayNoise, size_t nRow,
                                     float*  const arrayData,
                                     float         cutoff);
template
void array_remove_rezo<double>(const double* const arrayNoise, size_t nRow,
                                     double* const arrayData,
                                     double        cutoff);

} // namespace gem

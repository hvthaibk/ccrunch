/***********************************************************************
 *  File:       array_math.cpp
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
 * Math
 ****************************************/

// abs
template <typename T>
void array_math_abs(T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] = (T) std::abs(array[iRow]);
    }
}

template <>
void array_math_abs<uint8_t >(uint8_t*  const , size_t ) { }
template <>
void array_math_abs<uint16_t>(uint16_t* const , size_t ) { }
template <>
void array_math_abs<uint32_t>(uint32_t* const , size_t ) { }
template <>
void array_math_abs<uint64_t>(uint64_t* const , size_t ) { }

template <typename T>
void array_math_abs(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(arraySrc != arrayDst);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = (T) std::abs(arraySrc[iRow]);
    }
}

template <>
void array_math_abs<uint8_t >(uint8_t*  const , const uint8_t*  const , size_t ) { }
template <>
void array_math_abs<uint16_t>(uint16_t* const , const uint16_t* const , size_t ) { }
template <>
void array_math_abs<uint32_t>(uint32_t* const , const uint32_t* const , size_t ) { }
template <>
void array_math_abs<uint64_t>(uint64_t* const , const uint64_t* const , size_t ) { }

// instantiation
template
void array_math_abs<int8_t >(int8_t*  const array, size_t nRow);
template
void array_math_abs<int16_t>(int16_t* const array, size_t nRow);
template
void array_math_abs<int32_t>(int32_t* const array, size_t nRow);
template
void array_math_abs<int64_t>(int64_t* const array, size_t nRow);
template
void array_math_abs<float  >(float*   const array, size_t nRow);
template
void array_math_abs<double >(double*  const array, size_t nRow);
template
void array_math_abs<int8_t >(int8_t*  const arrayDst, const int8_t*  const arraySrc, size_t nRow);
template
void array_math_abs<int16_t>(int16_t* const arrayDst, const int16_t* const arraySrc, size_t nRow);
template
void array_math_abs<int32_t>(int32_t* const arrayDst, const int32_t* const arraySrc, size_t nRow);
template
void array_math_abs<int64_t>(int64_t* const arrayDst, const int64_t* const arraySrc, size_t nRow);
template
void array_math_abs<float  >(float*   const arrayDst, const float*   const arraySrc, size_t nRow);
template
void array_math_abs<double >(double*  const arrayDst, const double*  const arraySrc, size_t nRow);

// inverse
template <typename T>
void array_math_inv(T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] = -array[iRow];
    }
}

template <typename T>
void array_math_inv(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = -arraySrc[iRow];
    }
}

// instantiation
template
void array_math_inv<int32_t >(int32_t*  const array, size_t nRow);
template
void array_math_inv<uint32_t>(uint32_t* const array, size_t nRow);
template
void array_math_inv<int64_t >(int64_t*  const array, size_t nRow);
template
void array_math_inv<uint64_t>(uint64_t* const array, size_t nRow);
template
void array_math_inv<float   >(float*    const array, size_t nRow);
template
void array_math_inv<double  >(double*   const array, size_t nRow);
template
void array_math_inv<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void array_math_inv<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void array_math_inv<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void array_math_inv<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void array_math_inv<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void array_math_inv<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// norm (mean = 0, std = 1)
template <typename T>
void array_math_norm(T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    T         mean, std;

    mean = array_reduce_sum(array, nRow) / (T) nRow;

    array_math_sub(array, mean, nRow);

    std = (T) std::sqrt(array_reduce_sum2(array, nRow) / (T) nRow);

    array_math_div(array, std, nRow);
}

template <typename T>
void array_math_norm(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    T         mean, std;

    mean = array_reduce_sum(arraySrc, nRow) / (T) nRow;

    array_math_sub(arrayDst, arraySrc, mean, nRow);

    std = (T) std::sqrt(array_reduce_sum2(arrayDst, nRow) / (T) nRow);

    array_math_div(arrayDst, std, nRow);
}

// instantiation
template
void array_math_norm<int32_t >(int32_t*  const array, size_t nRow);
template
void array_math_norm<uint32_t>(uint32_t* const array, size_t nRow);
template
void array_math_norm<int64_t >(int64_t*  const array, size_t nRow);
template
void array_math_norm<uint64_t>(uint64_t* const array, size_t nRow);
template
void array_math_norm<float   >(float*    const array, size_t nRow);
template
void array_math_norm<double  >(double*   const array, size_t nRow);
template
void array_math_norm<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void array_math_norm<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void array_math_norm<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void array_math_norm<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void array_math_norm<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void array_math_norm<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// square
template <typename T>
void array_math_sqr(T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] *= array[iRow];
    }
}

template <typename T>
void array_math_sqr(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = arraySrc[iRow] * arraySrc[iRow];
    }
}

// instantiation
template
void array_math_sqr<int32_t >(int32_t*  const array, size_t nRow);
template
void array_math_sqr<uint32_t>(uint32_t* const array, size_t nRow);
template
void array_math_sqr<int64_t >(int64_t*  const array, size_t nRow);
template
void array_math_sqr<uint64_t>(uint64_t* const array, size_t nRow);
template
void array_math_sqr<float   >(float*    const array, size_t nRow);
template
void array_math_sqr<double  >(double*   const array, size_t nRow);
template
void array_math_sqr<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void array_math_sqr<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void array_math_sqr<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void array_math_sqr<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void array_math_sqr<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void array_math_sqr<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// square root
template <typename T>
void array_math_sqrt(T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] = (T) std::sqrt(array[iRow]);
    }
}

template <typename T>
void array_math_sqrt(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = (T) std::sqrt(arraySrc[iRow]);
    }
}

// instantiation
template
void array_math_sqrt<int32_t >(int32_t*  const array, size_t nRow);
template
void array_math_sqrt<uint32_t>(uint32_t* const array, size_t nRow);
template
void array_math_sqrt<int64_t >(int64_t*  const array, size_t nRow);
template
void array_math_sqrt<uint64_t>(uint64_t* const array, size_t nRow);
template
void array_math_sqrt<float   >(float*    const array, size_t nRow);
template
void array_math_sqrt<double  >(double*   const array, size_t nRow);
template
void array_math_sqrt<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void array_math_sqrt<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void array_math_sqrt<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void array_math_sqrt<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void array_math_sqrt<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void array_math_sqrt<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);

// add
template <typename T>
void array_math_add(T* const array, T value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] += value;
    }
}

template <typename T>
void array_math_add(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] += arraySrc[iRow];
    }
}

template <typename T>
void array_math_add(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = arraySrc[iRow] + value;
    }
}

template <typename T>
void array_math_add(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(arrayDst  != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = arraySrc1[iRow] + arraySrc2[iRow];
    }
}

// instantiation
template
void array_math_add<int32_t >(int32_t * const array, int32_t  value, size_t nRow);
template
void array_math_add<uint32_t>(uint32_t* const array, uint32_t value, size_t nRow);
template
void array_math_add<int64_t >(int64_t * const array, int64_t  value, size_t nRow);
template
void array_math_add<uint64_t>(uint64_t* const array, uint64_t value, size_t nRow);
template
void array_math_add<float   >(float*    const array, float    value, size_t nRow);
template
void array_math_add<double  >(double*   const array, double   value, size_t nRow);
template
void array_math_add<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void array_math_add<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void array_math_add<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void array_math_add<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void array_math_add<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void array_math_add<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);
template
void array_math_add<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, int32_t  value, size_t nRow);
template
void array_math_add<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, uint32_t value, size_t nRow);
template
void array_math_add<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, int64_t  value, size_t nRow);
template
void array_math_add<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, uint64_t value, size_t nRow);
template
void array_math_add<float   >(float*    const arrayDst, const float*    const arraySrc, float    value, size_t nRow);
template
void array_math_add<double  >(double*   const arrayDst, const double*   const arraySrc, double   value, size_t nRow);
template
void array_math_add<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow);
template
void array_math_add<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow);
template
void array_math_add<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow);
template
void array_math_add<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow);
template
void array_math_add<float   >(float*    const arrayDst, const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow);
template
void array_math_add<double  >(double*   const arrayDst, const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow);

// subtract
template <typename T>
void array_math_sub(T* const array, T value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] -= value;
    }
}

template <typename T>
void array_math_sub(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] -= arraySrc[iRow];
    }
}

template <typename T>
void array_math_sub(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = arraySrc[iRow] - value;
    }
}

template <typename T>
void array_math_sub(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(arrayDst  != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = arraySrc1[iRow] - arraySrc2[iRow];
    }
}

// instantiation
template
void array_math_sub<int32_t >(int32_t * const array, int32_t  value, size_t nRow);
template
void array_math_sub<uint32_t>(uint32_t* const array, uint32_t value, size_t nRow);
template
void array_math_sub<int64_t >(int64_t * const array, int64_t  value, size_t nRow);
template
void array_math_sub<uint64_t>(uint64_t* const array, uint64_t value, size_t nRow);
template
void array_math_sub<float   >(float*    const array, float    value, size_t nRow);
template
void array_math_sub<double  >(double*   const array, double   value, size_t nRow);
template
void array_math_sub<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void array_math_sub<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void array_math_sub<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void array_math_sub<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void array_math_sub<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void array_math_sub<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);
template
void array_math_sub<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, int32_t  value, size_t nRow);
template
void array_math_sub<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, uint32_t value, size_t nRow);
template
void array_math_sub<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, int64_t  value, size_t nRow);
template
void array_math_sub<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, uint64_t value, size_t nRow);
template
void array_math_sub<float   >(float*    const arrayDst, const float*    const arraySrc, float    value, size_t nRow);
template
void array_math_sub<double  >(double*   const arrayDst, const double*   const arraySrc, double   value, size_t nRow);
template
void array_math_sub<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow);
template
void array_math_sub<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow);
template
void array_math_sub<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow);
template
void array_math_sub<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow);
template
void array_math_sub<float   >(float*    const arrayDst, const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow);
template
void array_math_sub<double  >(double*   const arrayDst, const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow);

// multiply
template <typename T>
void array_math_mul(T* const array, T value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] *= value;
    }
}

template <typename T>
void array_math_mul(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] *= arraySrc[iRow];
    }
}

template <typename T>
void array_math_mul(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = arraySrc[iRow] * value;
    }
}

template <typename T>
void array_math_mul(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(arrayDst  != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = arraySrc1[iRow] * arraySrc2[iRow];
    }
}

// instantiation
template
void array_math_mul<int32_t >(int32_t * const array, int32_t  value, size_t nRow);
template
void array_math_mul<uint32_t>(uint32_t* const array, uint32_t value, size_t nRow);
template
void array_math_mul<int64_t >(int64_t * const array, int64_t  value, size_t nRow);
template
void array_math_mul<uint64_t>(uint64_t* const array, uint64_t value, size_t nRow);
template
void array_math_mul<float   >(float*    const array, float    value, size_t nRow);
template
void array_math_mul<double  >(double*   const array, double   value, size_t nRow);
template
void array_math_mul<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void array_math_mul<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void array_math_mul<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void array_math_mul<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void array_math_mul<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void array_math_mul<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);
template
void array_math_mul<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, int32_t  value, size_t nRow);
template
void array_math_mul<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, uint32_t value, size_t nRow);
template
void array_math_mul<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, int64_t  value, size_t nRow);
template
void array_math_mul<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, uint64_t value, size_t nRow);
template
void array_math_mul<float   >(float*    const arrayDst, const float*    const arraySrc, float    value, size_t nRow);
template
void array_math_mul<double  >(double*   const arrayDst, const double*   const arraySrc, double   value, size_t nRow);
template
void array_math_mul<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow);
template
void array_math_mul<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow);
template
void array_math_mul<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow);
template
void array_math_mul<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow);
template
void array_math_mul<float   >(float*    const arrayDst, const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow);
template
void array_math_mul<double  >(double*   const arrayDst, const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow);

// divide
template <typename T>
void array_math_div(T* const array, T value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        array[iRow] /= value;
    }
}

template <typename T>
void array_math_div(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] /= arraySrc[iRow];
    }
}

template <typename T>
void array_math_div(T* const arrayDst, const T* const arraySrc, T value, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = arraySrc[iRow] / value;
    }
}

template <typename T>
void array_math_div(T* const arrayDst, const T* const arraySrc1, const T* const arraySrc2, size_t nRow)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(arrayDst  != NULL);
    assert(nRow > 0);

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        arrayDst[iRow] = arraySrc1[iRow] / arraySrc2[iRow];
    }
}

// instantiation
template
void array_math_div<int32_t >(int32_t * const array, int32_t  value, size_t nRow);
template
void array_math_div<uint32_t>(uint32_t* const array, uint32_t value, size_t nRow);
template
void array_math_div<int64_t >(int64_t * const array, int64_t  value, size_t nRow);
template
void array_math_div<uint64_t>(uint64_t* const array, uint64_t value, size_t nRow);
template
void array_math_div<float   >(float*    const array, float    value, size_t nRow);
template
void array_math_div<double  >(double*   const array, double   value, size_t nRow);
template
void array_math_div<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void array_math_div<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void array_math_div<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void array_math_div<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void array_math_div<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void array_math_div<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);
template
void array_math_div<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, int32_t  value, size_t nRow);
template
void array_math_div<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, uint32_t value, size_t nRow);
template
void array_math_div<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, int64_t  value, size_t nRow);
template
void array_math_div<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, uint64_t value, size_t nRow);
template
void array_math_div<float   >(float*    const arrayDst, const float*    const arraySrc, float    value, size_t nRow);
template
void array_math_div<double  >(double*   const arrayDst, const double*   const arraySrc, double   value, size_t nRow);
template
void array_math_div<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow);
template
void array_math_div<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow);
template
void array_math_div<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow);
template
void array_math_div<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow);
template
void array_math_div<float   >(float*    const arrayDst, const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow);
template
void array_math_div<double  >(double*   const arrayDst, const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow);

} // namespace gem

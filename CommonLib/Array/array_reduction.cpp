/***********************************************************************
 *  File:       array_reduction.cpp
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
 * Reduction
 ****************************************/

// max
template <typename T>
T array_reduce_max(const T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    T         maxTmp = array[0];

    for (size_t iRow = 1; iRow < nRow; iRow++) {
        maxTmp = (array[iRow] > maxTmp) ? array[iRow] : maxTmp;
    }

    return maxTmp;
}

// instantiation
template
int32_t  array_reduce_max<int32_t >(const int32_t*  const array, size_t nRow);
template
uint32_t array_reduce_max<uint32_t>(const uint32_t* const array, size_t nRow);
template
int64_t  array_reduce_max<int64_t >(const int64_t*  const array, size_t nRow);
template
uint64_t array_reduce_max<uint64_t>(const uint64_t* const array, size_t nRow);
template
float    array_reduce_max<float   >(const float*    const array, size_t nRow);
template
double   array_reduce_max<double  >(const double*   const array, size_t nRow);

// min
template <typename T>
T array_reduce_min(const T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    T         mimTmp = array[0];

    for (size_t iRow = 1; iRow < nRow; iRow++) {
        mimTmp = (array[iRow] < mimTmp) ? array[iRow] : mimTmp;
    }

    return mimTmp;
}

// instantiation
template
int32_t  array_reduce_min<int32_t >(const int32_t*  const array, size_t nRow);
template
uint32_t array_reduce_min<uint32_t>(const uint32_t* const array, size_t nRow);
template
int64_t  array_reduce_min<int64_t >(const int64_t*  const array, size_t nRow);
template
uint64_t array_reduce_min<uint64_t>(const uint64_t* const array, size_t nRow);
template
float    array_reduce_min<float   >(const float*    const array, size_t nRow);
template
double   array_reduce_min<double  >(const double*   const array, size_t nRow);

// sum
template <typename T>
T array_reduce_sum(const T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    T        sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        sum += array[iRow];
    }

    return sum;
}

// instantiation
template
int32_t  array_reduce_sum<int32_t >(const int32_t*  const array, size_t nRow);
template
uint32_t array_reduce_sum<uint32_t>(const uint32_t* const array, size_t nRow);
template
int64_t  array_reduce_sum<int64_t >(const int64_t*  const array, size_t nRow);
template
uint64_t array_reduce_sum<uint64_t>(const uint64_t* const array, size_t nRow);
template
float    array_reduce_sum<float   >(const float*    const array, size_t nRow);
template
double   array_reduce_sum<double  >(const double*   const array, size_t nRow);

// sum of squares
template <typename T>
T array_reduce_sum2(const T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    T        sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        sum += array[iRow]*array[iRow];
    }

    return sum;
}

// instantiation
template
int32_t  array_reduce_sum2<int32_t >(const int32_t*  const array, size_t nRow);
template
uint32_t array_reduce_sum2<uint32_t>(const uint32_t* const array, size_t nRow);
template
int64_t  array_reduce_sum2<int64_t >(const int64_t*  const array, size_t nRow);
template
uint64_t array_reduce_sum2<uint64_t>(const uint64_t* const array, size_t nRow);
template
float    array_reduce_sum2<float   >(const float*    const array, size_t nRow);
template
double   array_reduce_sum2<double  >(const double*   const array, size_t nRow);

// mean
template <typename T>
T array_reduce_mean(const T* const array, size_t nRow)
{
    return array_reduce_sum(array, nRow) / (T) nRow;
}

// instantiation
template
int32_t  array_reduce_mean<int32_t >(const int32_t*  const array, size_t nRow);
template
uint32_t array_reduce_mean<uint32_t>(const uint32_t* const array, size_t nRow);
template
int64_t  array_reduce_mean<int64_t >(const int64_t*  const array, size_t nRow);
template
uint64_t array_reduce_mean<uint64_t>(const uint64_t* const array, size_t nRow);
template
float    array_reduce_mean<float   >(const float*    const array, size_t nRow);
template
double   array_reduce_mean<double  >(const double*   const array, size_t nRow);

// std
template <typename T>
T array_reduce_std(const T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    T        sum2 = 0;
    T        mean = array_reduce_mean(array, nRow);

    #pragma omp parallel for reduction(+:sum2)
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        sum2 += (array[iRow]-mean)*(array[iRow]-mean);
    }

    return (T) std::sqrt(sum2 / (T) nRow);
}

// instantiation
template
int32_t  array_reduce_std<int32_t >(const int32_t*  const array, size_t nRow);
template
uint32_t array_reduce_std<uint32_t>(const uint32_t* const array, size_t nRow);
template
int64_t  array_reduce_std<int64_t >(const int64_t*  const array, size_t nRow);
template
uint64_t array_reduce_std<uint64_t>(const uint64_t* const array, size_t nRow);
template
float    array_reduce_std<float   >(const float*    const array, size_t nRow);
template
double   array_reduce_std<double  >(const double*   const array, size_t nRow);

// index max
template <typename T>
size_t array_reduce_maxindx(const T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    T         maxTmp  = array[0];
    size_t    indxTmp = 0;

    for (size_t iRow = 1; iRow < nRow; iRow++) {
        if (array[iRow] > maxTmp) {
            maxTmp  = array[iRow];
            indxTmp = iRow;
        }
    }

    return indxTmp;
}

// instantiation
template
size_t array_reduce_maxindx<int32_t >(const int32_t*  const array, size_t nRow);
template
size_t array_reduce_maxindx<uint32_t>(const uint32_t* const array, size_t nRow);
template
size_t array_reduce_maxindx<int64_t >(const int64_t*  const array, size_t nRow);
template
size_t array_reduce_maxindx<uint64_t>(const uint64_t* const array, size_t nRow);
template
size_t array_reduce_maxindx<float   >(const float*    const array, size_t nRow);
template
size_t array_reduce_maxindx<double  >(const double*   const array, size_t nRow);

// index min
template <typename T>
size_t array_reduce_minindx(const T* const array, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    T         mimTmp  = array[0];
    size_t    indxTmp = 0;

    for (size_t iRow = 1; iRow < nRow; iRow++) {
        if (array[iRow] < mimTmp) {
            mimTmp  = array[iRow];
            indxTmp = iRow;
        }
    }

    return indxTmp;
}

// instantiation
template
size_t array_reduce_minindx<int32_t >(const int32_t*  const array, size_t nRow);
template
size_t array_reduce_minindx<uint32_t>(const uint32_t* const array, size_t nRow);
template
size_t array_reduce_minindx<int64_t >(const int64_t*  const array, size_t nRow);
template
size_t array_reduce_minindx<uint64_t>(const uint64_t* const array, size_t nRow);
template
size_t array_reduce_minindx<float   >(const float*    const array, size_t nRow);
template
size_t array_reduce_minindx<double  >(const double*   const array, size_t nRow);

// corr
template <typename T>
T array_reduce_corr(const T* const arraySrc1,
                    const T* const arraySrc2,
                    size_t nRow)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(nRow > 0);

    T        sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        sum += arraySrc1[iRow]*arraySrc2[iRow];
    }

    return sum;
}

// instantiation
template
int32_t  array_reduce_corr(const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow);
template
uint32_t array_reduce_corr(const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow);
template
int64_t  array_reduce_corr(const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow);
template
uint64_t array_reduce_corr(const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow);
template
float    array_reduce_corr(const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow);
template
double   array_reduce_corr(const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow);

// conv
template <typename T>
T array_reduce_conv(const T* const arraySrc1,
                    const T* const arraySrc2,
                    size_t nRow)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(nRow > 0);

    T        sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (size_t iRow = 0; iRow < nRow; iRow++) {
        sum += arraySrc1[iRow]*arraySrc2[nRow-1-iRow];
    }

    return sum;
}

// instantiation
template
int32_t  array_reduce_conv(const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow);
template
uint32_t array_reduce_conv(const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow);
template
int64_t  array_reduce_conv(const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow);
template
uint64_t array_reduce_conv(const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow);
template
float    array_reduce_conv(const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow);
template
double   array_reduce_conv(const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow);

// compare
template <typename T>
bool array_reduce_compare(const T* const arraySrc1, const T* const arraySrc2, size_t nRow, T threshold)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(nRow > 0);
    assert(std::abs(threshold) > 0);

    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (std::abs(arraySrc1[iRow] - arraySrc2[iRow]) >= std::abs(threshold)) {
            return false;
        }
    }

    return true;
}

template <>
bool array_reduce_compare<uint8_t>(const uint8_t* const arraySrc1, const uint8_t* const arraySrc2, size_t nRow, uint8_t threshold)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(nRow > 0);
    assert(std::abs((int32_t) threshold) > 0);

    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (std::abs((int32_t) arraySrc1[iRow] - (int32_t) arraySrc2[iRow]) >= std::abs((int32_t) threshold)) {
            return false;
        }
    }

    return true;
}

template <>
bool array_reduce_compare<uint16_t>(const uint16_t* const arraySrc1, const uint16_t* const arraySrc2, size_t nRow, uint16_t threshold)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(nRow > 0);
    assert(std::abs((int32_t) threshold) > 0);

    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (std::abs((int32_t) arraySrc1[iRow] - (int32_t) arraySrc2[iRow]) >= std::abs((int32_t) threshold)) {
            return false;
        }
    }

    return true;
}

template <>
bool array_reduce_compare<uint32_t>(const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow, uint32_t threshold)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(nRow > 0);
    assert(std::abs((int32_t) threshold) > 0);

    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (std::abs((int32_t) arraySrc1[iRow] - (int32_t) arraySrc2[iRow]) >= std::abs((int32_t) threshold)) {
            return false;
        }
    }

    return true;
}

template <>
bool array_reduce_compare<uint64_t>(const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow, uint64_t threshold)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
    assert(nRow > 0);
    assert(std::abs((int64_t) threshold) > 0);

    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (std::abs((int64_t) arraySrc1[iRow] - (int64_t) arraySrc2[iRow]) >= std::abs((int64_t) threshold)) {
            return false;
        }
    }

    return true;
}

// instantiation
template
bool array_reduce_compare<int8_t >(const int8_t*  const arraySrc1, const int8_t*  const arraySrc2, size_t nRow, int8_t  threshold);
template
bool array_reduce_compare<int16_t>(const int16_t* const arraySrc1, const int16_t* const arraySrc2, size_t nRow, int16_t threshold);
template
bool array_reduce_compare<int32_t>(const int32_t* const arraySrc1, const int32_t* const arraySrc2, size_t nRow, int32_t threshold);
template
bool array_reduce_compare<int64_t>(const int64_t* const arraySrc1, const int64_t* const arraySrc2, size_t nRow, int64_t threshold);
template
bool array_reduce_compare<float  >(const float*   const arraySrc1, const float*   const arraySrc2, size_t nRow, float   threshold);
template
bool array_reduce_compare<double >(const double*  const arraySrc1, const double*  const arraySrc2, size_t nRow, double  threshold);
template
bool array_reduce_compare<std::complex<float > >(const std::complex<float >* const arraySrc1,
                                                 const std::complex<float >* const arraySrc2,
                                                 size_t nRow, std::complex<float > threshold);
template
bool array_reduce_compare<std::complex<double> >(const std::complex<double>* const arraySrc1,
                                                 const std::complex<double>* const arraySrc2,
                                                 size_t nRow, std::complex<double> threshold);

} // namespace gem

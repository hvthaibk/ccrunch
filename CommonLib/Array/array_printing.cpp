/***********************************************************************
 *  File:       array_printing.cpp
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

// print 1D complex array
template <typename T>
void array_print(const std::complex<T>* const array,
                 size_t nRow,
                 int width,
                 int precision,
                 bool scientific)
{
    assert(array != NULL);
    assert(nRow > 0);

    size_t    i;

    std::cout.precision(precision);
    if (scientific) std::cout << std::scientific;
    else            std::cout << std::fixed;

    for (i = 0; i < nRow; i++) {
        std::cout << " {"
                  << std::setw(width) << std::real(array[i]) << ", "
                  << std::setw(width) << std::imag(array[i]) << "} " << std::endl;
    }

    std::cout << std::endl;
}

// instantiation
template
void array_print<float >(const std::complex<float >* const array, size_t nRow, int width, int precision, bool scientific);
template
void array_print<double>(const std::complex<double>* const array, size_t nRow, int width, int precision, bool scientific);

// print 1D real array
template <typename T>
void array_print(const T* const array,
                 size_t nRow,
                 int width,
                 int precision,
                 bool scientific )
{
    assert(array != NULL);
    assert(nRow > 0);

    size_t    i;

    std::cout.precision(precision);
    if (scientific) std::cout << std::scientific;
    else            std::cout << std::fixed;

    for (i = 0; i < nRow; i++) {
        std::cout << std::setw(width) << array[i];
    }

    std::cout << std::endl;
    std::cout << std::endl;
}

// instantiation
template
void array_print<uint8_t >(const uint8_t*  const array, size_t nRow, int width, int precision, bool scientific);
template
void array_print<int32_t >(const int32_t*  const array, size_t nRow, int width, int precision, bool scientific);
template
void array_print<uint32_t>(const uint32_t* const array, size_t nRow, int width, int precision, bool scientific);
template
void array_print<int64_t >(const int64_t*  const array, size_t nRow, int width, int precision, bool scientific);
template
void array_print<uint64_t>(const uint64_t* const array, size_t nRow, int width, int precision, bool scientific);
template
void array_print<float   >(const float*    const array, size_t nRow, int width, int precision, bool scientific);
template
void array_print<double  >(const double*   const array, size_t nRow, int width, int precision, bool scientific);

// print 2D complex array (*)
template <typename T>
void array_print(const std::complex<T>* const array,
                 size_t nRow, size_t nCol,
                 int width,
                 int precision,
                 bool scientific)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0);

    size_t    i, j, indx;

    std::cout.precision(precision);
    if (scientific) std::cout << std::scientific;
    else            std::cout << std::fixed;

    for (i = 0, indx = 0; i < nRow; i++) {
        for (j = 0; j < nCol; j++, indx++) {
            std::cout << " {"
                      << std::setw(width) << std::real(array[indx]) << ", "
                      << std::setw(width) << std::real(array[indx]) << "} ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

// instantiation
template
void array_print<float >(const std::complex<float >* const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<double>(const std::complex<double>* const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);

// print 2D complex array (**)
template <typename T>
void array_print(std::complex<T> const * const * const array,
                 size_t nRow, size_t nCol,
                 int width,
                 int precision,
                 bool scientific)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0);

    size_t    i, j;

    std::cout.precision(precision);
    if (scientific) std::cout << std::scientific;
    else            std::cout << std::fixed;

    for (i = 0; i < nRow; i++) {
        for (j = 0; j < nCol; j++) {
            std::cout << " {"
                      << std::setw(width) << std::real(array[i][j]) << ", "
                      << std::setw(width) << std::imag(array[i][j]) << "} ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

// instantiation
template
void array_print<float >(std::complex<float > const * const * const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<double>(std::complex<double> const * const * const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);

// print 2D real array (*)
template <typename T>
void array_print(const T* const array,
                 size_t nRow, size_t nCol,
                 int width,
                 int precision,
                 bool scientific)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0);

    size_t    i, j, indx;

    std::cout.precision(precision);
    if (scientific) std::cout << std::scientific;
    else            std::cout << std::fixed;

    for (i = 0, indx = 0; i < nRow; i++) {
        for (j = 0; j < nCol; j++, indx++) {
            std::cout << std::setw(width) << array[indx];
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

// instantiation
template
void array_print<uint8_t >(const uint8_t*  const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<int32_t >(const int32_t*  const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<uint32_t>(const uint32_t* const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<int64_t >(const int64_t*  const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<uint64_t>(const uint64_t* const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<float   >(const float*    const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<double  >(const double*   const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);

// print 2D real array (**)
template <typename T>
void array_print(T const * const * const array,
                 size_t nRow, size_t nCol,
                 int width,
                 int precision,
                 bool scientific)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0);

    size_t    i, j;

    std::cout.precision(precision);
    if (scientific) std::cout << std::scientific;
    else            std::cout << std::fixed;

    for (i = 0; i < nRow; i++) {
        for (j = 0; j < nCol; j++) {
            std::cout << std::setw(width) << array[i][j];
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

// instantiation
template
void array_print<uint8_t >(uint8_t  const * const * const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<int32_t >(int32_t  const * const * const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<uint32_t>(uint32_t const * const * const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<int64_t >(int64_t  const * const * const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<uint64_t>(uint64_t const * const * const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<float   >(float    const * const * const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);
template
void array_print<double  >(double   const * const * const array, size_t nRow, size_t nCol, int width, int precision, bool scientific);

// print 3D complex array (*)
template <typename T>
void array_print(const std::complex<T>* const array,
                 size_t nRow, size_t nCol, size_t nSec,
                 int width,
                 int precision,
                 bool scientific)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    size_t    i, j, k, indx;

    std::cout.precision(precision);
    if (scientific) std::cout << std::scientific;
    else            std::cout << std::fixed;

    for (k = 0; k < nSec; k++) {
        std::cout << "(:,:," << k << ")" << std::endl;
        for (i = 0; i < nRow; i++) {
            for (j = 0; j < nCol; j++) {
                indx = i*nCol*nSec + j*nSec + k;
                std::cout << " {"
                          << std::setw(width) << std::real(array[indx]) << ", "
                          << std::setw(width) << std::imag(array[indx]) << "} ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << std::endl;
}

// instantiation
template
void array_print<float >(const std::complex<float >* const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<double>(const std::complex<double>* const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);

// print 3D complex array (***)
template <typename T>
void array_print(std::complex<T> const * const * const * const array,
                 size_t nRow, size_t nCol, size_t nSec,
                 int width,
                 int precision,
                 bool scientific)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    size_t    i, j, k;

    std::cout.precision(precision);
    if (scientific) std::cout << std::scientific;
    else            std::cout << std::fixed;

    for (k = 0; k < nSec; k++) {
        std::cout << "(:,:," << k << ")" << std::endl;
        for (i = 0; i < nRow; i++) {
            for (j = 0; j < nCol; j++) {
                std::cout << " {"
                          << std::setw(width) << std::real(array[i][j][k]) << ", "
                          << std::setw(width) << std::imag(array[i][j][k]) << "} ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << std::endl;
}

// instantiation
template
void array_print<float >(std::complex<float > const * const * const * const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<double>(std::complex<double> const * const * const * const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);

// print 3D real array (*)
template <typename T>
void array_print(const T* const array,
                 size_t nRow, size_t nCol, size_t nSec,
                 int width,
                 int precision,
                 bool scientific)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    size_t    i, j, k, indx;

    std::cout.precision(precision);
    if (scientific) std::cout << std::scientific;
    else            std::cout << std::fixed;

    for (k = 0; k < nSec; k++) {
        std::cout << "(:,:," << k << ")" << std::endl;
        for (i = 0; i < nRow; i++) {
            for (j = 0; j < nCol; j++) {
                indx = i*nCol*nSec + j*nSec + k;
                std::cout << std::setw(width) << array[indx];
            }
            std::cout << std::endl;
        }
    }

    std::cout << std::endl;
}

// instantiation
template
void array_print<uint8_t >(const uint8_t*  const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<int32_t >(const int32_t*  const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<uint32_t>(const uint32_t* const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<int64_t >(const int64_t*  const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<uint64_t>(const uint64_t* const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<float   >(const float*    const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<double  >(const double*   const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);

// print 3D real array (***)
template <typename T>
void array_print(T const * const * const * const array,
                 size_t nRow, size_t nCol, size_t nSec,
                 int width,
                 int precision,
                 bool scientific)
{
    assert(array != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    size_t    i, j, k;

    std::cout.precision(precision);
    if (scientific) std::cout << std::scientific;
    else            std::cout << std::fixed;

    for (k = 0; k < nSec; k++) {
        std::cout << "(:,:," << k << ")" << std::endl;
        for (i = 0; i < nRow; i++) {
            for (j = 0; j < nCol; j++) {
                std::cout << std::setw(width) << array[i][j][k];
            }
            std::cout << std::endl;
        }
    }

    std::cout << std::endl;
}

// instantiation
template
void array_print<uint8_t >(uint8_t  const * const * const * const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<int32_t >(int32_t  const * const * const * const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<uint32_t>(uint32_t const * const * const * const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<int64_t >(int64_t  const * const * const * const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<uint64_t>(uint64_t const * const * const * const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<float   >(float    const * const * const * const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);
template
void array_print<double  >(double   const * const * const * const array, size_t nRow, size_t nCol, size_t nSec, int width, int precision, bool scientific);

} // namespace gem

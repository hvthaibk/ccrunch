/***********************************************************************
 *  File:       transform_expand_scaling.cpp
 *
 *  Purpose:    Implementation of expansion functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "transform.hpp"

namespace gem {

/**************************
 * Expansion
 *************************/

// 1D
template <typename T>
void transform_expand_scale_full(const T* const arraySrc, size_t  nRowSrc,
                                       T*&      arrayDst, size_t& nRowDst,
                                       T        factor)
{
    assert(arraySrc != NULL && arrayDst == NULL);
    assert(nRowSrc > 0);

    size_t    diagLength;
    size_t    nRowOff;

    diagLength = (size_t) std::ceil(factor*dist((T)0,(T)nRowSrc));

    nRowDst = ((diagLength-nRowSrc) % 2) ? diagLength+1 : diagLength;

    nRowDst = std::max(nRowDst, nRowSrc);

    nRowOff = (nRowDst-nRowSrc) / 2;

    array_new_zero(arrayDst, nRowDst);

    array_replace<T,T>(arrayDst, nRowDst,
                       arraySrc, nRowSrc,
                                 nRowOff);
}

// instantiation
template
void transform_expand_scale_full<float >(const float*  const arraySrc, size_t  nRowSrc,
                                               float*&       arrayDst, size_t& nRowDst,
                                               float         factor);
template
void transform_expand_scale_full<double>(const double* const arraySrc, size_t  nRowSrc,
                                               double*&      arrayDst, size_t& nRowDst,
                                               double        factor);

// 2D
template <typename T>
void transform_expand_scale_full(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                       T*&      arrayDst, size_t& nRowDst, size_t& nColDst,
                                       T        factor)
{
    assert(arraySrc != NULL && arrayDst == NULL);
    assert(nRowSrc > 0 && nColSrc > 0);

    size_t    diagLength;
    size_t    nRowOff, nColOff;

    diagLength = (size_t) std::ceil(factor*dist((T)0,(T)0,(T)nRowSrc,(T)nColSrc));

    nRowDst = ((diagLength-nRowSrc) % 2) ? diagLength+1 : diagLength;
    nColDst = ((diagLength-nColSrc) % 2) ? diagLength+1 : diagLength;

    nRowDst = std::max(nRowDst, nRowSrc);
    nColDst = std::max(nColDst, nColSrc);

    nRowOff = (nRowDst-nRowSrc) / 2;
    nColOff = (nColDst-nColSrc) / 2;

    array_new_zero(arrayDst, nRowDst*nColDst);

    array_replace<T,T>(arrayDst, nRowDst, nColDst,
                       arraySrc, nRowSrc, nColSrc,
                                 nRowOff, nColOff);
}

// instantiation
template
void transform_expand_scale_full<float >(const float*  const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                               float*&       arrayDst, size_t& nRowDst, size_t& nColDst,
                                               float         factor);
template
void transform_expand_scale_full<double>(const double* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                               double*&      arrayDst, size_t& nRowDst, size_t& nColDst,
                                               double        factor);

// 3D
template <typename T>
void transform_expand_scale_full(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                       T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst,
                                       T        factor)
{
    assert(arraySrc != NULL && arrayDst == NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);

    size_t    diagLength;
    size_t    nRowOff, nColOff, nSecOff;

    diagLength = (size_t) std::ceil(factor*dist((T)0,(T)0,(T)0,(T)nRowSrc,(T)nColSrc,(T)nSecSrc));

    nRowDst = ((diagLength-nRowSrc) % 2) ? diagLength+1 : diagLength;
    nColDst = ((diagLength-nColSrc) % 2) ? diagLength+1 : diagLength;
    nSecDst = ((diagLength-nSecSrc) % 2) ? diagLength+1 : diagLength;

    nRowDst = std::max(nRowDst, nRowSrc);
    nColDst = std::max(nColDst, nColSrc);
    nSecDst = std::max(nSecDst, nSecSrc);

    nRowOff = (nRowDst-nRowSrc) / 2;
    nColOff = (nColDst-nColSrc) / 2;
    nSecOff = (nSecDst-nSecSrc) / 2;

    array_new_zero(arrayDst, nRowDst*nColDst*nSecDst);

    array_replace<T,T>(arrayDst, nRowDst, nColDst, nSecDst,
                       arraySrc, nRowSrc, nColSrc, nSecSrc,
                                 nRowOff, nColOff, nSecOff);
}

// instantiation
template
void transform_expand_scale_full<float >(const float*  const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                               float*&       arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst,
                                               float         factor);
template
void transform_expand_scale_full<double>(const double* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                               double*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst,
                                               double        factor);

} // namespace gem

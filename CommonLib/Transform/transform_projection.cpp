/***********************************************************************
 *  File:       transform_projection.cpp
 *
 *  Purpose:    Implementation of transformation functions
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
 * Projection
 *************************/

// 2D
template <typename T>
void transform_project(const T* const arraySrc, size_t nRowSrc, size_t nColSrc,
                             T* const arrayDst,
                       const std::vector<T> arrayTheta,
                       eInter inter)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0);
    assert(nColSrc > 0);
    assert(arrayTheta.size() > 0);

    T*          arrayRot = NULL;

    array_new(arrayRot, nRowSrc*nColSrc);

    for (size_t i = 0; i < arrayTheta.size(); i++) {

        transform_rotate(arraySrc, arrayRot,
                         nRowSrc, nColSrc,
                         -arrayTheta[i],
                         inter);

        array_reducendim_sum(arrayRot, nRowSrc, nColSrc,
                             arrayDst+i*nColSrc,
                             REDUCE_NDIM_ROW);
    }

    array_delete(arrayRot);
}

// instantiation
template
void transform_project<float >(const float*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                     float*  const arrayDst,
                               const std::vector<float > arrayTheta,
                               eInter inter);
template
void transform_project<double>(const double* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                     double* const arrayDst,
                               const std::vector<double> arrayTheta,
                               eInter inter);

// 2D
template <typename T>
void transform_project(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                             T*&      arrayDst, size_t& nRowDst, size_t& nColDst,
                       const std::vector<T> arrayTheta,
                       eInter inter)
{
    assert(arraySrc != NULL);
    assert(arrayDst == NULL);
    assert(nRowSrc > 0);
    assert(nColSrc > 0);
    assert(arrayTheta.size() > 0);

    size_t      nRowExp, nColExp;
    T*          arrayExp = NULL;
    T*          arrayRot = NULL;

    transform_expand_rotate_full(arraySrc, nRowSrc, nColSrc,
                                 arrayExp, nRowExp, nColExp);
    array_new(arrayRot, nRowExp*nColExp);

    nRowDst = arrayTheta.size();
    nColDst = nColExp;
    array_new(arrayDst, nRowDst*nColDst);

    for (size_t i = 0; i < arrayTheta.size(); i++) {

        transform_rotate(arrayExp, arrayRot,
                         nRowExp, nColExp,
                         -arrayTheta[i],
                         inter);

        array_reducendim_sum(arrayRot, nRowExp, nColExp,
                             arrayDst+i*nColExp,
                             REDUCE_NDIM_ROW);
    }

    array_delete(arrayExp);
    array_delete(arrayRot);
}

// instantiation
template
void transform_project<float >(const float*  const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                     float*&       arrayDst, size_t& nRowDst, size_t& nColDst,
                               const std::vector<float > arrayTheta,
                               eInter inter);
template
void transform_project<double>(const double* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                     double*&      arrayDst, size_t& nRowDst, size_t& nColDst,
                               const std::vector<double> arrayTheta,
                               eInter inter);

// 3D
template <typename T>
void transform_project(const T* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                             T* const arrayDst,
                       const std::vector<T> arrayAlpha,
                       const std::vector<T> arrayBeta,
                       eInter inter)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRowSrc > 0);
    assert(nColSrc > 0);
    assert(nSecSrc > 0);
    assert(arrayAlpha.size() > 0);
    assert(arrayAlpha.size() == arrayBeta.size());

    T*          arrayRot = NULL;

    array_new(arrayRot, nRowSrc*nColSrc*nSecSrc);

    for (size_t i = 0; i < arrayAlpha.size(); i++) {

        transform_rotate(arraySrc, arrayRot,
                         nRowSrc, nColSrc, nSecSrc,
                         (T) 0, -arrayBeta[i], -arrayAlpha[i],
                         inter);

        array_reducendim_sum(arrayRot, nRowSrc, nColSrc, nSecSrc,
                             arrayDst+i*nColSrc*nSecSrc,
                             REDUCE_NDIM_ROW);
    }

    array_delete(arrayRot);
}

// instantiation
template
void transform_project<float >(const float*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                     float*  const arrayDst,
                               const std::vector<float > arrayAlpha,
                               const std::vector<float > arrayBeta,
                               eInter inter);
template
void transform_project<double>(const double* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                     double* const arrayDst,
                               const std::vector<double> arrayAlpha,
                               const std::vector<double> arrayBeta,
                               eInter inter);

// 3D
template <typename T>
void transform_project(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                             T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst,
                       const std::vector<T> arrayAlpha,
                       const std::vector<T> arrayBeta,
                       eInter inter)
{
    assert(arraySrc != NULL);
    assert(arrayDst == NULL);
    assert(nRowSrc > 0);
    assert(nColSrc > 0);
    assert(nSecSrc > 0);
    assert(arrayAlpha.size() > 0);
    assert(arrayAlpha.size() == arrayBeta.size());

    size_t      nRowExp, nColExp, nSecExp;
    T*          arrayExp = NULL;
    T*          arrayRot = NULL;

    transform_expand_rotate_full(arraySrc, nRowSrc, nColSrc, nSecSrc,
                                 arrayExp, nRowExp, nColExp, nSecExp);
    array_new(arrayRot, nRowExp*nColExp*nSecExp);

    nRowDst = arrayAlpha.size();
    nColDst = nColExp;
    nSecDst = nSecExp;
    array_new(arrayDst, nRowDst*nColDst*nSecDst);

    for (size_t i = 0; i < arrayAlpha.size(); i++) {

        transform_rotate(arrayExp, arrayRot,
                         nRowExp, nColExp, nSecExp,
                         (T) 0, -arrayBeta[i], -arrayAlpha[i],
                         inter);

        array_reducendim_sum(arrayRot, nRowExp, nColExp, nSecExp,
                             arrayDst+i*nColDst*nSecDst,
                             REDUCE_NDIM_ROW);
    }

    array_delete(arrayExp);
    array_delete(arrayRot);
}

// instantiation
template
void transform_project<float >(const float*  const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                     float*&       arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst,
                               const std::vector<float > arrayAlpha,
                               const std::vector<float > arrayBeta,
                               eInter inter);
template
void transform_project<double>(const double* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                     double*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst,
                               const std::vector<double> arrayAlpha,
                               const std::vector<double> arrayBeta,
                               eInter inter);

} // namespace gem

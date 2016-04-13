/***********************************************************************
 *  File:       transform_tracing.cpp
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
 * Tracing
 *************************/

// 2D
template <typename T>
void transform_trace(const T* const arraySrc, size_t nRowSrc, size_t nColSrc,
                           T* const arrayDst,
                     const std::vector<T> arrayTheta,
                     eTrace trace,
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

        switch (trace) {
            case TRACE_MAX:
                array_reducendim_max (arrayRot, nRowSrc, nColSrc,
                                      arrayDst+i*nColSrc,
                                      REDUCE_NDIM_ROW);
                break;
            case TRACE_MIN:
                array_reducendim_min (arrayRot, nRowSrc, nColSrc,
                                      arrayDst+i*nColSrc,
                                      REDUCE_NDIM_ROW);
                break;
            case TRACE_SUM:
                array_reducendim_sum (arrayRot, nRowSrc, nColSrc,
                                      arrayDst+i*nColSrc,
                                      REDUCE_NDIM_ROW);
                break;
            case TRACE_SUM2:
                array_reducendim_sum2(arrayRot, nRowSrc, nColSrc,
                                      arrayDst+i*nColSrc,
                                      REDUCE_NDIM_ROW);
                break;
            case TRACE_MEAN:
                array_reducendim_mean(arrayRot, nRowSrc, nColSrc,
                                      arrayDst+i*nColSrc,
                                      REDUCE_NDIM_ROW);
                break;
            case TRACE_STD:
                array_reducendim_std (arrayRot, nRowSrc, nColSrc,
                                      arrayDst+i*nColSrc,
                                      REDUCE_NDIM_ROW);
                break;
            default:
                ERROR("transform_trace", "unsupported tracing mode");
        }
    }

    array_delete(arrayRot);
}

// instantiation
template
void transform_trace<float >(const float*  const arraySrc, size_t nRowSrc, size_t nColSrc,
                                   float*  const arrayDst,
                             const std::vector<float > arrayTheta,
                             eTrace trace,
                             eInter inter);
template
void transform_trace<double>(const double* const arraySrc, size_t nRowSrc, size_t nColSrc,
                                   double* const arrayDst,
                             const std::vector<double> arrayTheta,
                             eTrace trace,
                             eInter inter);

// 3D
template <typename T>
void transform_trace(const T* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                           T* const arrayDst,
                     const std::vector<T> arrayAlpha,
                     const std::vector<T> arrayBeta,
                     const std::vector<T> arrayGamma,
                     eTrace trace,
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
                         arrayAlpha[i], arrayBeta[i], arrayGamma[i],
                         inter);

        switch (trace) {
            case TRACE_MAX:
                array_reducendim_max (arrayRot, nRowSrc, nColSrc, nSecSrc,
                                      arrayDst+i*nRowSrc*nColSrc,
                                      REDUCE_NDIM_SEC);
                break;
            case TRACE_MIN:
                array_reducendim_min (arrayRot, nRowSrc, nColSrc, nSecSrc,
                                      arrayDst+i*nRowSrc*nColSrc,
                                      REDUCE_NDIM_SEC);
                break;
            case TRACE_SUM:
                array_reducendim_sum (arrayRot, nRowSrc, nColSrc, nSecSrc,
                                      arrayDst+i*nRowSrc*nColSrc,
                                      REDUCE_NDIM_SEC);
                break;
            case TRACE_SUM2:
                array_reducendim_sum2(arrayRot, nRowSrc, nColSrc, nSecSrc,
                                      arrayDst+i*nRowSrc*nColSrc,
                                      REDUCE_NDIM_SEC);
                break;
            case TRACE_MEAN:
                array_reducendim_mean(arrayRot, nRowSrc, nColSrc, nSecSrc,
                                      arrayDst+i*nRowSrc*nColSrc,
                                      REDUCE_NDIM_SEC);
                break;
            case TRACE_STD:
                array_reducendim_std (arrayRot, nRowSrc, nColSrc, nSecSrc,
                                      arrayDst+i*nRowSrc*nColSrc,
                                      REDUCE_NDIM_SEC);
                break;
            default:
                ERROR("transform_trace", "unsupported tracing mode");
        }
    }

    array_delete(arrayRot);
}

// instantiation
template
void transform_trace<float >(const float*  const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                   float*  const arrayDst,
                             const std::vector<float > arrayAlpha,
                             const std::vector<float > arrayBeta,
                             const std::vector<float > arrayGamma,
                             eTrace trace,
                             eInter inter);
template
void transform_trace<double>(const double* const arraySrc, size_t nRowSrc, size_t nColSrc, size_t nSecSrc,
                                   double* const arrayDst,
                             const std::vector<double> arrayAlpha,
                             const std::vector<double> arrayBeta,
                             const std::vector<double> arrayGamma,
                             eTrace trace,
                             eInter inter);

} // namespace gem

/***********************************************************************
 *  File:       xcorr_direct.cpp
 *
 *  Purpose:    Implementation of xcorr-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "xcorr.hpp"

namespace gem {

// xcorr
template <typename T>
T xcorr_direct(const T* const tplData,
               const T* const refData,
               size_t length)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(length > 0);

    return array_reduce_corr(tplData, refData, length);
}

// instantiation
template
float  xcorr_direct<float >(const float*  const tplData,
                            const float*  const refData,
                            size_t length);
template
double xcorr_direct<double>(const double* const tplData,
                            const double* const refData,
                            size_t length);

// ecorr
template <typename T>
T ecorr_direct(const T* const tplData,
               const T* const refData,
               size_t length)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(length > 0);

    T     num, den;

    num = xcorr_direct(tplData, refData, length);

    den = std::sqrt(array_reduce_sum2(tplData, length) *
                    array_reduce_sum2(refData, length));

    return num / den;
}

// instantiation
template
float  ecorr_direct<float >(const float*  const tplData,
                            const float*  const refData,
                            size_t length);
template
double ecorr_direct<double>(const double* const tplData,
                            const double* const refData,
                            size_t length);

// ecorrm
template <typename T>
T ecorrm_direct(const T* const tplData,
                const T* const refData,
                const T* const mskData,
                size_t length)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(mskData != NULL);
    assert(length > 0);

    T     *tplTmp = NULL;
    T     *refTmp = NULL;
    T     corr;

    array_new(tplTmp, length);
    array_new(refTmp, length);

    array_memcpy(tplTmp, tplData, length);
    array_memcpy(refTmp, refData, length);

    xcorrNormTemplateECC(tplTmp, length, mskData);
    xcorrNormTemplateECC(refTmp, length, mskData);

    corr = xcorr_direct(tplTmp, refTmp, length);

    array_delete(tplTmp);
    array_delete(refTmp);

    return corr;
}

// instantiation
template
float  ecorrm_direct<float >(const float*  const tplData,
                             const float*  const refData,
                             const float*  const mskData,
                             size_t length);
template
double ecorrm_direct<double>(const double* const tplData,
                             const double* const refData,
                             const double* const mskData,
                             size_t length);

// normxcorr
template <typename T>
T normxcorr_direct(const T* const tplData,
                   const T* const refData,
                   size_t length)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(length > 0);

    T     *tplTmp = NULL;
    T     *refTmp = NULL;
    T     num, den;

    array_new(tplTmp, length);
    array_new(refTmp, length);

    array_math_sub(tplTmp, tplData, array_reduce_sum(tplData, length) / (T) length, length);
    array_math_sub(refTmp, refData, array_reduce_sum(refData, length) / (T) length, length);

    num = xcorr_direct(tplTmp, refTmp, length);

    den = std::sqrt(array_reduce_sum2(tplTmp, length) *
                    array_reduce_sum2(refTmp, length));

    array_delete(tplTmp);
    array_delete(refTmp);

    return num / den;
}

// instantiation
template
float  normxcorr_direct<float >(const float*  const tplData,
                                const float*  const refData,
                                size_t length);
template
double normxcorr_direct<double>(const double* const tplData,
                                const double* const refData,
                                size_t length);

// normxcorrm
template <typename T>
T normxcorrm_direct(const T* const tplData,
                    const T* const refData,
                    const T* const mskData,
                    size_t length)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(mskData != NULL);
    assert(length > 0);

    T     *tplTmp = NULL;
    T     *refTmp = NULL;
    T     mskArea, corr;

    array_new(tplTmp, length);
    array_new(refTmp, length);

    array_memcpy(tplTmp, tplData, length);
    array_memcpy(refTmp, refData, length);

    mskArea = (T) xcorrNormTemplateNCC(tplTmp, length, mskData);
    mskArea = (T) xcorrNormTemplateNCC(refTmp, length, mskData);

    corr = xcorr_direct(tplTmp, refTmp, length) / mskArea;

    array_delete(tplTmp);
    array_delete(refTmp);

    return corr;
}

// instantiation
template
float  normxcorrm_direct<float >(const float*  const tplData,
                                 const float*  const refData,
                                 const float*  const mskData,
                                 size_t length);
template
double normxcorrm_direct<double>(const double* const tplData,
                                 const double* const refData,
                                 const double* const mskData,
                                 size_t length);

// normxcorrmw
template <typename T>
T normxcorrmw_direct(const T* const tplData, const T* const tplWght,
                     const T* const refData, const T* const refWght,
                     const T* const mskData,
                     size_t length)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(mskData != NULL);
    assert(tplWght != NULL);
    assert(refWght != NULL);
    assert(length > 0);

    T     *tplTmp = NULL;
    T     *refTmp = NULL;
    T     *tplWghtTmp = NULL;
    T     *refWghtTmp = NULL;
    T     corr;

    array_new(tplTmp, length);
    array_new(refTmp, length);
    array_new(tplWghtTmp, length);
    array_new(refWghtTmp, length);

    array_memcpy(tplTmp, tplData, length);
    array_memcpy(refTmp, refData, length);
    array_memcpy(tplWghtTmp, tplWght, length);
    array_memcpy(refWghtTmp, refWght, length);

    xcorrNormTemplateWCC1(tplTmp, length, mskData, tplWghtTmp);
    xcorrNormTemplateWCC1(refTmp, length, mskData, refWghtTmp);

    corr = xcorr_direct(tplTmp, refTmp, length);

    array_delete(tplTmp);
    array_delete(refTmp);
    array_delete(tplWghtTmp);
    array_delete(refWghtTmp);

    return corr;
}

// instantiation
template
float  normxcorrmw_direct<float >(const float*  const tplData, const float*  const tplWght,
                                  const float*  const refData, const float*  const refWght,
                                  const float*  const mskData,
                                  size_t length);
template
double normxcorrmw_direct<double>(const double* const tplData, const double* const tplWght,
                                  const double* const refData, const double* const refWght,
                                  const double* const mskData,
                                  size_t length);

} // namespace gem

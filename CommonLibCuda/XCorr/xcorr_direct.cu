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

#include "xcorr.cuh"

namespace gem {

// xcorr
template <typename T>
T cuda_xcorr_direct(const T* const tplData,
                    const T* const refData,
                    size_t length)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(length > 0);

    return cuda_array_reduce_corr(tplData, refData, length);
}

// instantiation
template
float  cuda_xcorr_direct<float >(const float*  const tplData,
                                 const float*  const refData,
                                 size_t length);
template
double cuda_xcorr_direct<double>(const double* const tplData,
                                 const double* const refData,
                                 size_t length);

// ecorr
template <typename T>
T cuda_ecorr_direct(const T* const tplData,
                    const T* const refData,
                    size_t length)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(length > 0);

    T     num, den;

    num = cuda_xcorr_direct(tplData, refData, length);

    den = std::sqrt(cuda_array_reduce_sum2(tplData, length) *
                    cuda_array_reduce_sum2(refData, length));

    return num / den;
}

// instantiation
template
float  cuda_ecorr_direct<float >(const float*  const tplData,
                                 const float*  const refData,
                                 size_t length);
template
double cuda_ecorr_direct<double>(const double* const tplData,
                                 const double* const refData,
                                 size_t length);

// ecorrm
template <typename T>
T cuda_ecorrm_direct(const T* const tplData,
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

    cuda_arrayDev_new(tplTmp, length);
    cuda_arrayDev_new(refTmp, length);

    cuda_array_memcpy_d2d(tplTmp, tplData, length);
    cuda_array_memcpy_d2d(refTmp, refData, length);

    cuda_xcorrNormTemplateECC(tplTmp, length, mskData);
    cuda_xcorrNormTemplateECC(refTmp, length, mskData);

    corr = cuda_xcorr_direct(tplTmp, refTmp, length);

    cuda_arrayDev_delete(tplTmp);
    cuda_arrayDev_delete(refTmp);

    return corr;
}

// instantiation
template
float  cuda_ecorrm_direct<float >(const float*  const tplData,
                                  const float*  const refData,
                                  const float*  const mskData,
                                  size_t length);
template
double cuda_ecorrm_direct<double>(const double* const tplData,
                                  const double* const refData,
                                  const double* const mskData,
                                  size_t length);

// normxcorr
template <typename T>
T cuda_normxcorr_direct(const T* const tplData,
                        const T* const refData,
                        size_t length)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(length > 0);

    T     *tplTmp = NULL;
    T     *refTmp = NULL;
    T     num, den;

    cuda_arrayDev_new(tplTmp, length);
    cuda_arrayDev_new(refTmp, length);

    cuda_array_math_sub(tplTmp, tplData, cuda_array_reduce_sum(tplData, length) / (T) length, length);
    cuda_array_math_sub(refTmp, refData, cuda_array_reduce_sum(refData, length) / (T) length, length);

    num = cuda_xcorr_direct(tplTmp, refTmp, length);

    den = std::sqrt(cuda_array_reduce_sum2(tplTmp, length) *
                    cuda_array_reduce_sum2(refTmp, length));

    cuda_arrayDev_delete(tplTmp);
    cuda_arrayDev_delete(refTmp);

    return num / den;
}

// instantiation
template
float  cuda_normxcorr_direct<float >(const float*  const tplData,
                                     const float*  const refData,
                                     size_t length);
template
double cuda_normxcorr_direct<double>(const double* const tplData,
                                     const double* const refData,
                                     size_t length);

// normxcorrm
template <typename T>
T cuda_normxcorrm_direct(const T* const tplData,
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

    cuda_arrayDev_new(tplTmp, length);
    cuda_arrayDev_new(refTmp, length);

    cuda_array_memcpy_d2d(tplTmp, tplData, length);
    cuda_array_memcpy_d2d(refTmp, refData, length);

    mskArea = (T) cuda_xcorrNormTemplateNCC(tplTmp, length, mskData);
    mskArea = (T) cuda_xcorrNormTemplateNCC(refTmp, length, mskData);

    corr = cuda_xcorr_direct(tplTmp, refTmp, length) / mskArea;

    cuda_arrayDev_delete(tplTmp);
    cuda_arrayDev_delete(refTmp);

    return corr;
}

// instantiation
template
float  cuda_normxcorrm_direct<float >(const float*  const tplData,
                                      const float*  const refData,
                                      const float*  const mskData,
                                      size_t length);
template
double cuda_normxcorrm_direct<double>(const double* const tplData,
                                      const double* const refData,
                                      const double* const mskData,
                                      size_t length);

// normxcorrmw
template <typename T>
T cuda_normxcorrmw_direct(const T* const tplData, const T* const tplWght,
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

    cuda_arrayDev_new(tplTmp, length);
    cuda_arrayDev_new(refTmp, length);
    cuda_arrayDev_new(tplWghtTmp, length);
    cuda_arrayDev_new(refWghtTmp, length);

    cuda_array_memcpy_d2d(tplTmp, tplData, length);
    cuda_array_memcpy_d2d(refTmp, refData, length);
    cuda_array_memcpy_d2d(tplWghtTmp, tplWght, length);
    cuda_array_memcpy_d2d(refWghtTmp, refWght, length);

    cuda_xcorrNormTemplateWCC1(tplTmp, length, mskData, tplWghtTmp);
    cuda_xcorrNormTemplateWCC1(refTmp, length, mskData, refWghtTmp);

    corr = cuda_xcorr_direct(tplTmp, refTmp, length);

    cuda_arrayDev_delete(tplTmp);
    cuda_arrayDev_delete(refTmp);
    cuda_arrayDev_delete(tplWghtTmp);
    cuda_arrayDev_delete(refWghtTmp);

    return corr;
}

// instantiation
template
float  cuda_normxcorrmw_direct<float >(const float*  const tplData, const float*  const tplWght,
                                       const float*  const refData, const float*  const refWght,
                                       const float*  const mskData,
                                       size_t length);
template
double cuda_normxcorrmw_direct<double>(const double* const tplData, const double* const tplWght,
                                       const double* const refData, const double* const refWght,
                                       const double* const mskData,
                                       size_t length);

} // namespace gem

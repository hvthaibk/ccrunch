/***********************************************************************
 *  File:       xcorr_normxcorrmw.cu
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

/***********************************************************************
 * Normalized cross-correlation with mask + weight
 **********************************************************************/

// 1D
template <typename T>
void cuda_normxcorrmw(T *dev_tplData, size_t tplNrow, T *dev_tplWght,
                      T *dev_refData, size_t refNrow, T *dev_refWght,
                      T *dev_resData,
                      T *dev_mskData,
                      eXCorrRes shape)
{
    assert(dev_tplData != NULL && dev_refData != NULL && dev_resData != NULL);
    assert(dev_tplWght != NULL && dev_refWght != NULL && dev_mskData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    T         *dev_tplDataTmp = NULL,
              *dev_refDataTmp = NULL,
              *dev_tplWghtTmp = NULL,
              *dev_refWghtTmp = NULL,
              *dev_fftTmp1 = NULL,
              *dev_fftTmp2 = NULL;
    size_t    resNrow = 0, tplSize, refSize, resSize;

    if (shape != XCORR_RES_FULL && shape != XCORR_RES_VALID) {
        shape = XCORR_RES_FULL;
    }
    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow - tplNrow + 1;
            break;
        default:
            ERROR("cuda_normxcorrmw", "unsupported xcorr mode");
    }

    tplSize = tplNrow;
    refSize = refNrow;
    resSize = resNrow;

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_arrayDev_new(dev_refDataTmp, refSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp, dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataTmp, dev_refData, refSize);

    cuda_arrayDev_new(dev_tplWghtTmp, tplSize);
    cuda_arrayDev_new(dev_refWghtTmp, tplSize);
    cuda_array_memcpy_d2d(dev_tplWghtTmp, dev_tplWght, tplSize);
    cuda_array_memcpy_d2d(dev_refWghtTmp, dev_refWght, tplSize);

    // weight and template normalization
    cuda_xcorrNormTemplateWCC2(dev_tplDataTmp, tplSize,
                               dev_mskData,
                               dev_tplWghtTmp,
                               dev_refWghtTmp);

    // allocate temporary arrays
    cuda_arrayDev_new(dev_fftTmp1, resSize);
    cuda_arrayDev_new(dev_fftTmp2, resSize);

    // numerator
    cuda_xcorr(dev_tplDataTmp, tplNrow,
               dev_refDataTmp, refNrow,
               dev_resData,
               shape);

    // denominator
    cuda_xcorr(dev_refWghtTmp, tplNrow,
               dev_refDataTmp, refNrow,
               dev_fftTmp1,
               shape);

    cuda_array_math_sqr(dev_refDataTmp, refSize);
    cuda_xcorr(dev_refWghtTmp, tplNrow,
               dev_refDataTmp, refNrow,
               dev_fftTmp2,
               shape);

    // NCC
    cuda_xcorrCombineResultWCC(dev_resData,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               cuda_array_reduce_sum(dev_tplDataTmp, tplSize),
                               resSize);

    // deallocate
    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_refDataTmp);
    cuda_arrayDev_delete(dev_tplWghtTmp);
    cuda_arrayDev_delete(dev_refWghtTmp);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// instantiation
template
void cuda_normxcorrmw<float >(float  *dev_tplData, size_t tplNrow, float *dev_tplWght,
                              float  *dev_refData, size_t refNrow, float *dev_refWght,
                              float  *dev_resData,
                              float  *dev_mskData,
                              eXCorrRes shape);
template
void cuda_normxcorrmw<double>(double *dev_tplData, size_t tplNrow, double *dev_tplWght,
                              double *dev_refData, size_t refNrow, double *dev_refWght,
                              double *dev_resData,
                              double *dev_mskData,
                              eXCorrRes shape);

// 2D
template <typename T>
void cuda_normxcorrmw(T *dev_tplData, size_t tplNrow, size_t tplNcol, T *dev_tplWght,
                      T *dev_refData, size_t refNrow, size_t refNcol, T *dev_refWght,
                      T *dev_resData,
                      T *dev_mskData,
                      eXCorrRes shape)
{
    assert(dev_tplData != NULL && dev_refData != NULL && dev_resData != NULL);
    assert(dev_tplWght != NULL && dev_refWght != NULL && dev_mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    T         *dev_tplDataTmp = NULL,
              *dev_refDataTmp = NULL,
              *dev_tplWghtTmp = NULL,
              *dev_refWghtTmp = NULL,
              *dev_fftTmp1 = NULL,
              *dev_fftTmp2 = NULL;
    size_t    resNrow = 0, resNcol = 0, tplSize, refSize, resSize;

    if (shape != XCORR_RES_FULL && shape != XCORR_RES_VALID) {
        shape = XCORR_RES_FULL;
    }
    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            resNcol = refNcol + tplNcol - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow - tplNrow + 1;
            resNcol = refNcol - tplNcol + 1;
            break;
        default:
            ERROR("cuda_normxcorrmw", "unsupported xcorr mode");
    }

    tplSize = tplNrow * tplNcol;
    refSize = refNrow * refNcol;
    resSize = resNrow * resNcol;

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_arrayDev_new(dev_refDataTmp, refSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp, dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataTmp, dev_refData, refSize);

    cuda_arrayDev_new(dev_tplWghtTmp, tplSize);
    cuda_arrayDev_new(dev_refWghtTmp, tplSize);
    cuda_array_memcpy_d2d(dev_tplWghtTmp, dev_tplWght, tplSize);
    cuda_array_memcpy_d2d(dev_refWghtTmp, dev_refWght, tplSize);

    // weight and template normalization
    cuda_xcorrNormTemplateWCC2(dev_tplDataTmp, tplSize,
                               dev_mskData,
                               dev_tplWghtTmp,
                               dev_refWghtTmp);

    // allocate temporary arrays
    cuda_arrayDev_new(dev_fftTmp1, resSize);
    cuda_arrayDev_new(dev_fftTmp2, resSize);

    // numerator
    cuda_xcorr(dev_tplDataTmp, tplNrow, tplNcol,
               dev_refDataTmp, refNrow, refNcol,
               dev_resData,
               shape);

    // denominator
    cuda_xcorr(dev_refWghtTmp, tplNrow, tplNcol,
               dev_refDataTmp, refNrow, refNcol,
               dev_fftTmp1,
               shape);

    cuda_array_math_sqr(dev_refDataTmp, refSize);
    cuda_xcorr(dev_refWghtTmp, tplNrow, tplNcol,
               dev_refDataTmp, refNrow, refNcol,
               dev_fftTmp2,
               shape);

    // NCC
    cuda_xcorrCombineResultWCC(dev_resData,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               cuda_array_reduce_sum(dev_tplDataTmp, tplSize),
                               resSize);

    // deallocate
    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_refDataTmp);
    cuda_arrayDev_delete(dev_tplWghtTmp);
    cuda_arrayDev_delete(dev_refWghtTmp);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// instantiation
template
void cuda_normxcorrmw<float >(float  *dev_tplData, size_t tplNrow, size_t tplNcol, float *dev_tplWght,
                              float  *dev_refData, size_t refNrow, size_t refNcol, float *dev_refWght,
                              float  *dev_resData,
                              float  *dev_mskData,
                              eXCorrRes shape);
template
void cuda_normxcorrmw<double>(double *dev_tplData, size_t tplNrow, size_t tplNcol, double *dev_tplWght,
                              double *dev_refData, size_t refNrow, size_t refNcol, double *dev_refWght,
                              double *dev_resData,
                              double *dev_mskData,
                              eXCorrRes shape);

// 3D
template <typename T>
void cuda_normxcorrmw(T *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec, T *dev_tplWght,
                      T *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec, T *dev_refWght,
                      T *dev_resData,
                      T *dev_mskData,
                      eXCorrRes shape)
{
    assert(dev_tplData != NULL && dev_refData != NULL && dev_resData != NULL);
    assert(dev_tplWght != NULL && dev_refWght != NULL && dev_mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    T         *dev_tplDataTmp = NULL,
              *dev_refDataTmp = NULL,
              *dev_tplWghtTmp = NULL,
              *dev_refWghtTmp = NULL,
              *dev_fftTmp1 = NULL,
              *dev_fftTmp2 = NULL;
    size_t    resNrow = 0, resNcol = 0, resNsec = 0, tplSize, refSize, resSize;

    if (shape != XCORR_RES_FULL && shape != XCORR_RES_VALID) {
        shape = XCORR_RES_FULL;
    }
    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            resNcol = refNcol + tplNcol - 1;
            resNsec = refNsec + tplNsec - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow - tplNrow + 1;
            resNcol = refNcol - tplNcol + 1;
            resNsec = refNsec - tplNsec + 1;
            break;
        default:
            ERROR("cuda_normxcorrmw", "unsupported xcorr mode");
    }

    tplSize = tplNrow * tplNcol * tplNsec;
    refSize = refNrow * refNcol * refNsec;
    resSize = resNrow * resNcol * resNsec;

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_arrayDev_new(dev_refDataTmp, refSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp, dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataTmp, dev_refData, refSize);

    cuda_arrayDev_new(dev_tplWghtTmp, tplSize);
    cuda_arrayDev_new(dev_refWghtTmp, tplSize);
    cuda_array_memcpy_d2d(dev_tplWghtTmp, dev_tplWght, tplSize);
    cuda_array_memcpy_d2d(dev_refWghtTmp, dev_refWght, tplSize);

    // weight and template normalization
    cuda_xcorrNormTemplateWCC2(dev_tplDataTmp, tplSize,
                               dev_mskData,
                               dev_tplWghtTmp,
                               dev_refWghtTmp);

    // allocate temporary arrays
    cuda_arrayDev_new(dev_fftTmp1, resSize);
    cuda_arrayDev_new(dev_fftTmp2, resSize);

    // numerator
    cuda_xcorr(dev_tplDataTmp, tplNrow, tplNcol, tplNsec,
               dev_refDataTmp, refNrow, refNcol, refNsec,
               dev_resData,
               shape);

    // denominator
    cuda_xcorr(dev_refWghtTmp, tplNrow, tplNcol, tplNsec,
               dev_refDataTmp, refNrow, refNcol, refNsec,
               dev_fftTmp1,
               shape);

    cuda_array_math_sqr(dev_refDataTmp, refSize);
    cuda_xcorr(dev_refWghtTmp, tplNrow, tplNcol, tplNsec,
               dev_refDataTmp, refNrow, refNcol, refNsec,
               dev_fftTmp2,
               shape);

    // NCC
    cuda_xcorrCombineResultWCC(dev_resData,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               cuda_array_reduce_sum(dev_tplDataTmp, tplSize),
                               resSize);

    // deallocate
    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_refDataTmp);
    cuda_arrayDev_delete(dev_tplWghtTmp);
    cuda_arrayDev_delete(dev_refWghtTmp);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// instantiation
template
void cuda_normxcorrmw<float >(float  *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec, float *dev_tplWght,
                              float  *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec, float *dev_refWght,
                              float  *dev_resData,
                              float  *dev_mskData,
                              eXCorrRes shape);
template
void cuda_normxcorrmw<double>(double *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec, double *dev_tplWght,
                              double *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec, double *dev_refWght,
                              double *dev_resData,
                              double *dev_mskData,
                              eXCorrRes shape);

} // namespace gem

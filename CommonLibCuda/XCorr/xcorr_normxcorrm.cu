/***********************************************************************
 *  File:       xcorr_normxcorrm.cu
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
 * Normalized cross-correlation with mask
 **********************************************************************/

// 1D
template <typename T>
void cuda_normxcorrm(T *dev_tplData, size_t tplNrow,
                     T *dev_refData, size_t refNrow,
                     T *dev_resData,
                     T *dev_mskData,
                     eXCorrRes shape)
{
    assert(dev_tplData != NULL);
    assert(dev_refData != NULL);
    assert(dev_resData != NULL);
    assert(dev_mskData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    T         *dev_tplDataTmp = NULL,
              *dev_refDataTmp = NULL,
              *dev_fftTmp1 = NULL,
              *dev_fftTmp2 = NULL;
    size_t    resNrow, tplSize, refSize, resSize;
    size_t    mskArea;

    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow - tplNrow + 1;
            break;
        default:
            ERROR("cuda_normxcorrm", "unsupported xcorr mode");
    }

    tplSize = tplNrow;
    refSize = refNrow;
    resSize = resNrow;

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_arrayDev_new(dev_refDataTmp, refSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp, dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataTmp, dev_refData, refSize);

    // mask size and template normalization
    mskArea = cuda_xcorrNormTemplateNCC(dev_tplDataTmp, tplSize, dev_mskData);

    // allocate temporary arrays
    cuda_arrayDev_new(dev_fftTmp1, resSize);
    cuda_arrayDev_new(dev_fftTmp2, resSize);

    // numerator
    cuda_xcorr(dev_tplDataTmp, tplNrow,
               dev_refDataTmp, refNrow,
               dev_resData,
               shape);

    // denominator
    cuda_xcorr(dev_mskData,    tplNrow,
               dev_refDataTmp, refNrow,
               dev_fftTmp1,
               shape);

    cuda_array_math_sqr(dev_refDataTmp, refSize);
    cuda_xcorr(dev_mskData,    tplNrow,
               dev_refDataTmp, refNrow,
               dev_fftTmp2,
               shape);

    // NCC
    cuda_xcorrCombineResultNCC(dev_resData,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               mskArea,
                               resSize);

    // deallocate
    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_refDataTmp);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// instantiation
template
void cuda_normxcorrm<float >(float  *dev_tplData, size_t tplNrow,
                             float  *dev_refData, size_t refNrow,
                             float  *dev_resData,
                             float  *dev_mskData,
                             eXCorrRes shape);
template
void cuda_normxcorrm<double>(double *dev_tplData, size_t tplNrow,
                             double *dev_refData, size_t refNrow,
                             double *dev_resData,
                             double *dev_mskData,
                             eXCorrRes shape);

// 2D
template <typename T>
void cuda_normxcorrm(T *dev_tplData, size_t tplNrow, size_t tplNcol,
                     T *dev_refData, size_t refNrow, size_t refNcol,
                     T *dev_resData,
                     T *dev_mskData,
                     eXCorrRes shape)
{
    assert(dev_tplData != NULL);
    assert(dev_refData != NULL);
    assert(dev_resData != NULL);
    assert(dev_mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    T         *dev_tplDataTmp = NULL,
              *dev_refDataTmp = NULL,
              *dev_fftTmp1 = NULL,
              *dev_fftTmp2 = NULL;
    size_t    resNrow, resNcol, tplSize, refSize, resSize;
    size_t    mskArea;

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
            ERROR("cuda_normxcorrm", "unsupported xcorr mode");
    }

    tplSize = tplNrow * tplNcol;
    refSize = refNrow * refNcol;
    resSize = resNrow * resNcol;

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_arrayDev_new(dev_refDataTmp, refSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp, dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataTmp, dev_refData, refSize);

    // mask size and template normalization
    mskArea = cuda_xcorrNormTemplateNCC(dev_tplDataTmp, tplSize, dev_mskData);

    // allocate temporary arrays
    cuda_arrayDev_new(dev_fftTmp1, resSize);
    cuda_arrayDev_new(dev_fftTmp2, resSize);

    // numerator
    cuda_xcorr(dev_tplDataTmp, tplNrow, tplNcol,
               dev_refDataTmp, refNrow, refNcol,
               dev_resData,
               shape);

    // denominator
    cuda_xcorr(dev_mskData,    tplNrow, tplNcol,
               dev_refDataTmp, refNrow, refNcol,
               dev_fftTmp1,
               shape);

    cuda_array_math_sqr(dev_refDataTmp, refSize);
    cuda_xcorr(dev_mskData,    tplNrow, tplNcol,
               dev_refDataTmp, refNrow, refNcol,
               dev_fftTmp2,
               shape);

    // NCC
    cuda_xcorrCombineResultNCC(dev_resData,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               mskArea,
                               resSize);

    // deallocate
    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_refDataTmp);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// instantiation
template
void cuda_normxcorrm<float >(float  *dev_tplData, size_t tplNrow, size_t tplNcol,
                             float  *dev_refData, size_t refNrow, size_t refNcol,
                             float  *dev_resData,
                             float  *dev_mskData,
                             eXCorrRes shape);
template
void cuda_normxcorrm<double>(double *dev_tplData, size_t tplNrow, size_t tplNcol,
                             double *dev_refData, size_t refNrow, size_t refNcol,
                             double *dev_resData,
                             double *dev_mskData,
                             eXCorrRes shape);

// 3D
template <typename T>
void cuda_normxcorrm(T *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                     T *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                     T *dev_resData,
                     T *dev_mskData,
                     eXCorrRes shape)
{
    assert(dev_tplData != NULL);
    assert(dev_refData != NULL);
    assert(dev_resData != NULL);
    assert(dev_mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    T         *dev_tplDataTmp = NULL,
              *dev_refDataTmp = NULL,
              *dev_fftTmp1 = NULL,
              *dev_fftTmp2 = NULL;
    size_t    resNrow, resNcol, resNsec, tplSize, refSize, resSize;
    size_t    mskArea;

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
            ERROR("cuda_normxcorrm", "unsupported xcorr mode");
    }

    tplSize = tplNrow * tplNcol * tplNsec;
    refSize = refNrow * refNcol * refNsec;
    resSize = resNrow * resNcol * resNsec;

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_arrayDev_new(dev_refDataTmp, refSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp, dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataTmp, dev_refData, refSize);

    // mask size and template normalization
    mskArea = cuda_xcorrNormTemplateNCC(dev_tplDataTmp, tplSize, dev_mskData);

    // allocate temporary arrays
    cuda_arrayDev_new(dev_fftTmp1, resSize);
    cuda_arrayDev_new(dev_fftTmp2, resSize);

    // numerator
    cuda_xcorr(dev_tplDataTmp, tplNrow, tplNcol, tplNsec,
               dev_refDataTmp, refNrow, refNcol, refNsec,
               dev_resData,
               shape);

    // denominator
    cuda_xcorr(dev_mskData,    tplNrow, tplNcol, tplNsec,
               dev_refDataTmp, refNrow, refNcol, refNsec,
               dev_fftTmp1,
               shape);

    cuda_array_math_sqr(dev_refDataTmp, refSize);
    cuda_xcorr(dev_mskData,    tplNrow, tplNcol, tplNsec,
               dev_refDataTmp, refNrow, refNcol, refNsec,
               dev_fftTmp2,
               shape);

    // NCC
    cuda_xcorrCombineResultNCC(dev_resData,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               mskArea,
                               resSize);

    // deallocate
    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_refDataTmp);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// instantiation
template
void cuda_normxcorrm<float >(float  *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                             float  *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                             float  *dev_resData,
                             float  *dev_mskData,
                             eXCorrRes shape);
template
void cuda_normxcorrm<double>(double *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                             double *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                             double *dev_resData,
                             double *dev_mskData,
                             eXCorrRes shape);

} // namespace gem

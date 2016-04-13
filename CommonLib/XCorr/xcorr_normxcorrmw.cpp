/***********************************************************************
 *  File:       xcorr_normxcorrmw.cpp
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

/***************************************************
 * Normalized cross-correlation with mask + weight
 **************************************************/

// 1D
template <typename T>
void normxcorrmw(T *tplData, size_t tplNrow, T *tplWght,
                 T *refData, size_t refNrow, T *refWght,
                 T *resData,
                 T *mskData,
                 eXCorrRes shape)
{
    assert(tplData != NULL && refData != NULL && resData != NULL);
    assert(tplWght != NULL && refWght != NULL && mskData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    T         *tplDataTmp = NULL,
              *refDataTmp = NULL,
              *tplWghtTmp = NULL,
              *refWghtTmp = NULL,
              *fftTmp1 = NULL,
              *fftTmp2 = NULL;
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
            ERROR("normxcorrmw", "unsupported xcorr mode");
    }

    tplSize = tplNrow;
    refSize = refNrow;
    resSize = resNrow;

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    array_new(refDataTmp, refSize);
    memcpy(tplDataTmp, tplData, tplSize*sizeof(T));
    memcpy(refDataTmp, refData, refSize*sizeof(T));

    array_new(tplWghtTmp, tplSize);
    array_new(refWghtTmp, tplSize);
    memcpy(tplWghtTmp, tplWght, tplSize*sizeof(T));
    memcpy(refWghtTmp, refWght, tplSize*sizeof(T));

    // weight and template normalization
    xcorrNormTemplateWCC2(tplDataTmp, tplSize,
                          mskData,
                          tplWghtTmp,
                          refWghtTmp);

    // allocate temporary arrays
    array_new(fftTmp1, resSize);
    array_new(fftTmp2, resSize);

    // numerator
    xcorr(tplDataTmp, tplNrow,
          refDataTmp, refNrow,
          resData,
          shape);

    // denominator
    xcorr(refWghtTmp, tplNrow,
          refDataTmp, refNrow,
          fftTmp1,
          shape);

    array_math_sqr(refDataTmp, refSize);
    xcorr(refWghtTmp, tplNrow,
          refDataTmp, refNrow,
          fftTmp2,
          shape);

    // NCC
    xcorrCombineResultWCC(resData,
                          fftTmp1,
                          fftTmp2,
                          array_reduce_sum(tplDataTmp, tplSize),
                          resSize);

    // deallocate
    array_delete(tplDataTmp);
    array_delete(refDataTmp);
    array_delete(tplWghtTmp);
    array_delete(refWghtTmp);
    array_delete(fftTmp1);
    array_delete(fftTmp2);
}

// instantiation
template
void normxcorrmw<float >(float  *tplData, size_t tplNrow, float  *tplWght,
                         float  *refData, size_t refNrow, float  *refWght,
                         float  *resData,
                         float  *mskData,
                         eXCorrRes shape);
template
void normxcorrmw<double>(double *tplData, size_t tplNrow, double *tplWght,
                         double *refData, size_t refNrow, double *refWght,
                         double *resData,
                         double *mskData,
                         eXCorrRes shape);

// 2D
template <typename T>
void normxcorrmw(T *tplData, size_t tplNrow, size_t tplNcol, T *tplWght,
                 T *refData, size_t refNrow, size_t refNcol, T *refWght,
                 T *resData,
                 T *mskData,
                 eXCorrRes shape)
{
    assert(tplData != NULL && refData != NULL && resData != NULL);
    assert(tplWght != NULL && refWght != NULL && mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    T         *tplDataTmp = NULL,
              *refDataTmp = NULL,
              *tplWghtTmp = NULL,
              *refWghtTmp = NULL,
              *fftTmp1 = NULL,
              *fftTmp2 = NULL;
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
            ERROR("normxcorrmw", "unsupported xcorr mode");
    }

    tplSize = tplNrow * tplNcol;
    refSize = refNrow * refNcol;
    resSize = resNrow * resNcol;

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    array_new(refDataTmp, refSize);
    memcpy(tplDataTmp, tplData, tplSize*sizeof(T));
    memcpy(refDataTmp, refData, refSize*sizeof(T));

    array_new(tplWghtTmp, tplSize);
    array_new(refWghtTmp, tplSize);
    memcpy(tplWghtTmp, tplWght, tplSize*sizeof(T));
    memcpy(refWghtTmp, refWght, tplSize*sizeof(T));

    // weight and template normalization
    xcorrNormTemplateWCC2(tplDataTmp, tplSize,
                          mskData,
                          tplWghtTmp,
                          refWghtTmp);

    // allocate temporary arrays
    array_new(fftTmp1, resSize);
    array_new(fftTmp2, resSize);

    // numerator
    xcorr(tplDataTmp, tplNrow, tplNcol,
          refDataTmp, refNrow, refNcol,
          resData,
          shape);

    // denominator
    xcorr(refWghtTmp, tplNrow, tplNcol,
          refDataTmp, refNrow, refNcol,
          fftTmp1,
          shape);

    array_math_sqr(refDataTmp, refSize);
    xcorr(refWghtTmp, tplNrow, tplNcol,
          refDataTmp, refNrow, refNcol,
          fftTmp2,
          shape);

    // NCC
    xcorrCombineResultWCC(resData,
                          fftTmp1,
                          fftTmp2,
                          array_reduce_sum(tplDataTmp, tplSize),
                          resSize);

    // deallocate
    array_delete(tplDataTmp);
    array_delete(refDataTmp);
    array_delete(tplWghtTmp);
    array_delete(refWghtTmp);
    array_delete(fftTmp1);
    array_delete(fftTmp2);
}

// instantiation
template
void normxcorrmw<float >(float  *tplData, size_t tplNrow, size_t tplNcol, float  *tplWght,
                         float  *refData, size_t refNrow, size_t refNcol, float  *refWght,
                         float  *resData,
                         float  *mskData,
                         eXCorrRes shape);
template
void normxcorrmw<double>(double *tplData, size_t tplNrow, size_t tplNcol, double *tplWght,
                         double *refData, size_t refNrow, size_t refNcol, double *refWght,
                         double *resData,
                         double *mskData,
                         eXCorrRes shape);

// 3D
template <typename T>
void normxcorrmw(T *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec, T *tplWght,
                 T *refData, size_t refNrow, size_t refNcol, size_t refNsec, T *refWght,
                 T *resData,
                 T *mskData,
                 eXCorrRes shape)
{
    assert(tplData != NULL && refData != NULL && resData != NULL);
    assert(tplWght != NULL && refWght != NULL && mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    T         *tplDataTmp = NULL,
              *refDataTmp = NULL,
              *tplWghtTmp = NULL,
              *refWghtTmp = NULL,
              *fftTmp1 = NULL,
              *fftTmp2 = NULL;
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
            ERROR("normxcorrmw", "unsupported xcorr mode");
    }

    tplSize = tplNrow * tplNcol * tplNsec;
    refSize = refNrow * refNcol * refNsec;
    resSize = resNrow * resNcol * resNsec;

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    array_new(refDataTmp, refSize);
    memcpy(tplDataTmp, tplData, tplSize*sizeof(T));
    memcpy(refDataTmp, refData, refSize*sizeof(T));

    array_new(tplWghtTmp, tplSize);
    array_new(refWghtTmp, tplSize);
    memcpy(tplWghtTmp, tplWght, tplSize*sizeof(T));
    memcpy(refWghtTmp, refWght, tplSize*sizeof(T));

    // weight and template normalization
    xcorrNormTemplateWCC2(tplDataTmp, tplSize,
                          mskData,
                          tplWghtTmp,
                          refWghtTmp);

    // allocate temporary arrays
    array_new(fftTmp1, resSize);
    array_new(fftTmp2, resSize);

    // numerator
    xcorr(tplDataTmp, tplNrow, tplNcol, tplNsec,
          refDataTmp, refNrow, refNcol, refNsec,
          resData,
          shape);

    // denominator
    xcorr(refWghtTmp, tplNrow, tplNcol, tplNsec,
          refDataTmp, refNrow, refNcol, refNsec,
          fftTmp1,
          shape);

    array_math_sqr(refDataTmp, refSize);
    xcorr(refWghtTmp, tplNrow, tplNcol, tplNsec,
          refDataTmp, refNrow, refNcol, refNsec,
          fftTmp2,
          shape);

    // NCC
    xcorrCombineResultWCC(resData,
                          fftTmp1,
                          fftTmp2,
                          array_reduce_sum(tplDataTmp, tplSize),
                          resSize);

    // deallocate
    array_delete(tplDataTmp);
    array_delete(refDataTmp);
    array_delete(tplWghtTmp);
    array_delete(refWghtTmp);
    array_delete(fftTmp1);
    array_delete(fftTmp2);
}

// instantiation
template
void normxcorrmw<float >(float  *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec, float  *tplWght,
                         float  *refData, size_t refNrow, size_t refNcol, size_t refNsec, float  *refWght,
                         float  *resData,
                         float  *mskData,
                         eXCorrRes shape);
template
void normxcorrmw<double>(double *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec, double *tplWght,
                         double *refData, size_t refNrow, size_t refNcol, size_t refNsec, double *refWght,
                         double *resData,
                         double *mskData,
                         eXCorrRes shape);

} // namespace gem

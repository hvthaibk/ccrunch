/***********************************************************************
 *  File:       xcorr_regnormxcorrm.cpp
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
 * Masked registration using NCC
 **************************************************/

// 1D
template <typename T>
void regnormxcorrm(const T* const tplData, size_t tplNrow,
                   const T* const refData, size_t refNrow,
                   const T* const tplMask,
                   const T* const refMask,
                         T* const resCorr,
                         T* const resOlap,
                         T  overlapRatio)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(tplMask != NULL);
    assert(refMask != NULL);
    assert(resCorr != NULL);
    assert(resOlap != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);
    assert(overlapRatio > 0 && overlapRatio <= 1);

    T         *tplData2  = NULL,
              *refData2  = NULL;
    T         *numTpl    = NULL, *numRef = NULL;
    T         *denTpl    = NULL, *denRef = NULL;
    T         *resCutoff = NULL;
    size_t    resNrow;
    size_t    tplSize, refSize, resSize;

    resNrow = refNrow + tplNrow - 1;

    tplSize = tplNrow;
    refSize = refNrow;
    resSize = resNrow;

    // allocate
    array_new(tplData2,  tplSize);
    array_new(refData2,  refSize);
    array_new(numTpl,    resSize);
    array_new(numRef,    resSize);
    array_new(denTpl,    resSize);
    array_new(denRef,    resSize);
    array_new(resCutoff, resSize);

    // overlap matrix
    xcorr(tplMask, tplNrow,
          refMask, refNrow,
          resOlap,
          XCORR_RES_FULL);

    // numerator
    xcorr(tplData, tplNrow,
          refData, refNrow,
          resCorr,
          XCORR_RES_FULL);

    xcorr(tplData, tplNrow,
          refMask, refNrow,
          numTpl,
          XCORR_RES_FULL);

    xcorr(tplMask, tplNrow,
          refData, refNrow,
          numRef,
          XCORR_RES_FULL);

    // denumerator
    array_math_sqr(tplData2, tplData, tplSize);
    array_math_sqr(refData2, refData, refSize);

    xcorr(tplData2, tplNrow,
          refMask,  refNrow,
          denTpl,
          XCORR_RES_FULL);

    xcorr(tplMask,  tplNrow,
          refData2, refNrow,
          denRef,
          XCORR_RES_FULL);

    // NCC
    xcorrCombineResultNReg(resCorr, resOlap, numTpl, numRef, denTpl, denRef, resSize, overlapRatio);

    // deallocate
    array_delete(tplData2);
    array_delete(refData2);
    array_delete(numTpl);
    array_delete(numRef);
    array_delete(denTpl);
    array_delete(denRef);
    array_delete(resCutoff);
}

// instantiation
template
void regnormxcorrm<float >(const float*  const tplData, size_t tplNrow,
                           const float*  const refData, size_t refNrow,
                           const float*  const tplMask,
                           const float*  const refMask,
                                 float*  const resCorr,
                                 float*  const resOlap,
                                 float   overlapRatio);
template
void regnormxcorrm<double>(const double* const tplData, size_t tplNrow,
                           const double* const refData, size_t refNrow,
                           const double* const tplMask,
                           const double* const refMask,
                                 double* const resCorr,
                                 double* const resOlap,
                                 double  overlapRatio);

// 2D
template <typename T>
void regnormxcorrm(const T* const tplData, size_t tplNrow, size_t tplNcol,
                   const T* const refData, size_t refNrow, size_t refNcol,
                   const T* const tplMask,
                   const T* const refMask,
                         T* const resCorr,
                         T* const resOlap,
                         T  overlapRatio)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(tplMask != NULL);
    assert(refMask != NULL);
    assert(resCorr != NULL);
    assert(resOlap != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);
    assert(overlapRatio > 0 && overlapRatio <= 1);

    T         *tplData2  = NULL,
              *refData2  = NULL;
    T         *numTpl    = NULL, *numRef = NULL;
    T         *denTpl    = NULL, *denRef = NULL;
    T         *resCutoff = NULL;
    size_t    resNrow, resNcol;
    size_t    tplSize, refSize, resSize;

    resNrow = refNrow + tplNrow - 1;
    resNcol = refNcol + tplNcol - 1;

    tplSize = tplNrow * tplNcol;
    refSize = refNrow * refNcol;
    resSize = resNrow * resNcol;

    // allocate
    array_new(tplData2,  tplSize);
    array_new(refData2,  refSize);
    array_new(numTpl,    resSize);
    array_new(numRef,    resSize);
    array_new(denTpl,    resSize);
    array_new(denRef,    resSize);
    array_new(resCutoff, resSize);

    // overlap matrix
    xcorr(tplMask, tplNrow, tplNcol,
          refMask, refNrow, refNcol,
          resOlap,
          XCORR_RES_FULL);

    // numerator
    xcorr(tplData, tplNrow, tplNcol,
          refData, refNrow, refNcol,
          resCorr,
          XCORR_RES_FULL);

    xcorr(tplData, tplNrow, tplNcol,
          refMask, refNrow, refNcol,
          numTpl,
          XCORR_RES_FULL);

    xcorr(tplMask, tplNrow, tplNcol,
          refData, refNrow, refNcol,
          numRef,
          XCORR_RES_FULL);

    // denumerator
    array_math_sqr(tplData2, tplData, tplSize);
    array_math_sqr(refData2, refData, refSize);

    xcorr(tplData2, tplNrow, tplNcol,
          refMask,  refNrow, refNcol,
          denTpl,
          XCORR_RES_FULL);

    xcorr(tplMask,  tplNrow, tplNcol,
          refData2, refNrow, refNcol,
          denRef,
          XCORR_RES_FULL);

    // NCC
    xcorrCombineResultNReg(resCorr, resOlap, numTpl, numRef, denTpl, denRef, resSize, overlapRatio);

    // deallocate
    array_delete(tplData2);
    array_delete(refData2);
    array_delete(numTpl);
    array_delete(numRef);
    array_delete(denTpl);
    array_delete(denRef);
    array_delete(resCutoff);
}

// instantiation
template
void regnormxcorrm<float >(const float*  const tplData, size_t tplNrow, size_t tplNcol,
                           const float*  const refData, size_t refNrow, size_t refNcol,
                           const float*  const tplMask,
                           const float*  const refMask,
                                 float*  const resCorr,
                                 float*  const resOlap,
                                 float   overlapRatio);
template
void regnormxcorrm<double>(const double* const tplData, size_t tplNrow, size_t tplNcol,
                           const double* const refData, size_t refNrow, size_t refNcol,
                           const double* const tplMask,
                           const double* const refMask,
                                 double* const resCorr,
                                 double* const resOlap,
                                 double  overlapRatio);

// 3D
template <typename T>
void regnormxcorrm(const T* const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                   const T* const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                   const T* const tplMask,
                   const T* const refMask,
                         T* const resCorr,
                         T* const resOlap,
                         T  overlapRatio)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(tplMask != NULL);
    assert(refMask != NULL);
    assert(resCorr != NULL);
    assert(resOlap != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);
    assert(overlapRatio > 0 && overlapRatio <= 1);

    T         *tplData2  = NULL,
              *refData2  = NULL;
    T         *numTpl    = NULL, *numRef = NULL;
    T         *denTpl    = NULL, *denRef = NULL;
    T         *resCutoff = NULL;
    size_t    resNrow, resNcol, resNsec;
    size_t    tplSize, refSize, resSize;

    resNrow = refNrow + tplNrow - 1;
    resNcol = refNcol + tplNcol - 1;
    resNsec = refNsec + tplNsec - 1;

    tplSize = tplNrow * tplNcol * tplNsec;
    refSize = refNrow * refNcol * refNsec;
    resSize = resNrow * resNcol * resNsec;

    // allocate
    array_new(tplData2,  tplSize);
    array_new(refData2,  refSize);
    array_new(numTpl,    resSize);
    array_new(numRef,    resSize);
    array_new(denTpl,    resSize);
    array_new(denRef,    resSize);
    array_new(resCutoff, resSize);

    // overlap matrix
    xcorr(tplMask, tplNrow, tplNcol, tplNsec,
          refMask, refNrow, refNcol, refNsec,
          resOlap,
          XCORR_RES_FULL);

    // numerator
    xcorr(tplData, tplNrow, tplNcol, tplNsec,
          refData, refNrow, refNcol, refNsec,
          resCorr,
          XCORR_RES_FULL);

    xcorr(tplData, tplNrow, tplNcol, tplNsec,
          refMask, refNrow, refNcol, refNsec,
          numTpl,
          XCORR_RES_FULL);

    xcorr(tplMask, tplNrow, tplNcol, tplNsec,
          refData, refNrow, refNcol, refNsec,
          numRef,
          XCORR_RES_FULL);

    // denumerator
    array_math_sqr(tplData2, tplData, tplSize);
    array_math_sqr(refData2, refData, refSize);

    xcorr(tplData2, tplNrow, tplNcol, tplNsec,
          refMask,  refNrow, refNcol, refNsec,
          denTpl,
          XCORR_RES_FULL);

    xcorr(tplMask,  tplNrow, tplNcol, tplNsec,
          refData2, refNrow, refNcol, refNsec,
          denRef,
          XCORR_RES_FULL);

    // NCC
    xcorrCombineResultNReg(resCorr, resOlap, numTpl, numRef, denTpl, denRef, resSize, overlapRatio);

    // deallocate
    array_delete(tplData2);
    array_delete(refData2);
    array_delete(numTpl);
    array_delete(numRef);
    array_delete(denTpl);
    array_delete(denRef);
    array_delete(resCutoff);
}

// instantiation
template
void regnormxcorrm<float >(const float*  const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                           const float*  const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                           const float*  const tplMask,
                           const float*  const refMask,
                                 float*  const resCorr,
                                 float*  const resOlap,
                                 float   overlapRatio);
template
void regnormxcorrm<double>(const double* const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                           const double* const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                           const double* const tplMask,
                           const double* const refMask,
                                 double* const resCorr,
                                 double* const resOlap,
                                 double  overlapRatio);

} // namespace gem

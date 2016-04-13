/***********************************************************************
 *  File:       xcorr.hpp
 *
 *  Purpose:    Header file for xcorr-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_LIB_XCORR_HPP__
#define __GEM_LIB_XCORR_HPP__

#include "array.hpp"

namespace gem {

enum eXCorrRes
{
    XCORR_RES_FULL  = 0,
    XCORR_RES_SAME  = 1,
    XCORR_RES_VALID = 2
};

enum eXCorrMerge
{
    XCORR_MERGE_NEGATIVE = 0,       //     accept negative peaks
    XCORR_MERGE_POSITIVE = 1        // not accept negative peaks
};

enum eCCmode
{
    CC_MODE_XCORR = 0,      //                   cross-correlation
    CC_MODE_ECORR = 1,      // energy normalized cross-correlation
    CC_MODE_NCORR = 2,      //        normalized cross-correlation
    CC_MODE_WCORR = 3,      // weight normalized cross-correlation
    CC_MODE_SCORR = 4,      // surface correlation
    CC_MODE_CCORR = 5,      // combination of correaltions
    CC_MODE_MINFO = 6       // mutual information
};

/***************************************************
 * Helper functions for computing correlation
 **************************************************/

// modulation
template <typename T1, typename T2>
void pcorrModulateAndNormalize(const T1* const tplFFT,
                               const T1* const refFFT,
                                     T1* const resFFT,
                               size_t fftSize,
                               T2     nElements);

template <typename T1, typename T2>
void xcorrModulateAndNormalize(const T1* const tplFFT,
                               const T1* const refFFT,
                                     T1* const resFFT,
                               size_t fftSize,
                               T2     nElements);

// template normalization
template <typename T>
void xcorrNormTemplateECC(T* const tplData, size_t tplSize);

template <typename T>
void xcorrNormTemplateECC(      T* const tplData, size_t tplSize,
                          const T* const mskData);

template <typename T>
size_t xcorrNormTemplateNCC(T* const tplData, size_t tplSize);

template <typename T>
size_t xcorrNormTemplateNCC(      T* const tplData, size_t tplSize,
                            const T* const mskData);

template <typename T>
void xcorrNormTemplateWCC1(T* const tplData, size_t tplSize,
                           T* const tplWght);

template <typename T>
void xcorrNormTemplateWCC2(T* const tplData, size_t tplSize,
                           T* const tplWght,
                           T* const refWght);

template <typename T>
void xcorrNormTemplateWCC1(      T* const tplData, size_t tplSize,
                           const T* const mskData,
                                 T* const tplWght);

template <typename T>
void xcorrNormTemplateWCC2(      T* const tplData, size_t tplSize,
                           const T* const mskData,
                                 T* const tplWght,
                                 T* const refWght);

// combine correlations
template <typename T>
void xcorrCombineResultCCC(      T* const resData,
                           const T* const corrData1, T weight1,
                           const T* const corrData2, T weight2,
                           const T* const corrData3, T weight3,
                           size_t resSize);

template <typename T>
void xcorrCombineResultECC(      T* const resData,
                           const T* const fftTmp1,
                           size_t resSize);

template <typename T>
void xcorrCombineResultNCC(      T* const resData,
                           const T* const fftTmp1,
                           const T* const fftTmp2,
                           size_t mskArea,
                           size_t resSize);

template <typename T>
void xcorrCombineResultWCC(      T* const resData,
                           const T* const fftTmp1,
                           const T* const fftTmp2,
                           T      tplSum,
                           size_t resSize);

template <typename T>
void xcorrCombineResultNReg(      T* const resCorr,
                            const T* const resOlap,
                            const T* const numTpl,
                            const T* const numRef,
                            const T* const denTpl,
                            const T* const denRef,
                            size_t resSize,
                            T      overlapRatio);

// merge result
template <typename T>
void xcorrMergeResult(const T*      const resData,
                            T*      const resDataAbsMax,
                            size_t* const resDataMaxInd,
                      size_t        resSize,
                      size_t        indx,
                      eXCorrMerge   bAbs);

template <typename T>
void xcorrMergeResultGlobal(const T*      const resDataAbsMax,
                            const size_t* const resDataMaxInd,
                                  T*      const resDataAbsMaxGlobal,
                                  size_t* const resDataMaxIndGlobal,
                            size_t        resSize,
                            eXCorrMerge   bAbs);

/***************************************************
 * Phase-correlation
 **************************************************/

// size(tplData) = size(refData) = size(resData) = refNrow
void pcorr(const float*  const tplData,
           const float*  const refData,
                 float*  const resData,
           size_t refNrow);

void pcorr(const double* const tplData,
           const double* const refData,
                 double* const resData,
           size_t refNrow);

// size(tplData) = size(refData) = size(resData) = [refNrow refNcol]
void pcorr(const float*  const tplData,
           const float*  const refData,
                 float*  const resData,
           size_t refNrow, size_t refNcol);

void pcorr(const double* const tplData,
           const double* const refData,
                 double* const resData,
           size_t refNrow, size_t refNcol);

// size(tplData) = size(refData) = size(resData) = [refNrow refNcol refNsec]
void pcorr(const float*  const tplData,
           const float*  const refData,
                 float*  const resData,
           size_t refNrow, size_t refNcol, size_t refNsec);

void pcorr(const double* const tplData,
           const double* const refData,
                 double* const resData,
           size_t refNrow, size_t refNcol, size_t refNsec);

/***************************************************
 * Direct correlation
 **************************************************/

// xcorr
template <typename T>
T xcorr_direct(const T* const tplData,
               const T* const refData,
               size_t length);

// ecorr
template <typename T>
T ecorr_direct(const T* const tplData,
               const T* const refData,
               size_t length);

// ecorrm
template <typename T>
T ecorrm_direct(const T* const tplData,
                const T* const refData,
                const T* const mskData,
                size_t length);

// normxcorr
template <typename T>
T normxcorr_direct(const T* const tplData,
                   const T* const refData,
                   size_t length);

// normxcorrm
template <typename T>
T normxcorrm_direct(const T* const tplData,
                    const T* const refData,
                    const T* const mskData,
                    size_t length);

// normxcorrmw
template <typename T>
T normxcorrmw_direct(const T* const tplData, const T* const tplWght,
                     const T* const refData, const T* const refWght,
                     const T* const mskData,
                     size_t length);

/***************************************************
 * Cross-correlation
 **************************************************/

// 1D float
void xcorr(const float*  const tplData, size_t tplNrow,
           const float*  const refData, size_t refNrow,
                 float*  const resData,
           eXCorrRes shape = XCORR_RES_FULL);

// 1D double
void xcorr(const double* const tplData, size_t tplNrow,
           const double* const refData, size_t refNrow,
                 double* const resData,
           eXCorrRes shape = XCORR_RES_FULL);

// 2D float
void xcorr(const float*  const tplData, size_t tplNrow, size_t tplNcol,
           const float*  const refData, size_t refNrow, size_t refNcol,
                 float*  const resData,
           eXCorrRes shape = XCORR_RES_FULL);

// 2D double
void xcorr(const double* const tplData, size_t tplNrow, size_t tplNcol,
           const double* const refData, size_t refNrow, size_t refNcol,
                 double* const resData,
           eXCorrRes shape = XCORR_RES_FULL);

// 3D float
void xcorr(const float*  const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
           const float*  const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                 float*  const resData,
           eXCorrRes shape = XCORR_RES_FULL);

// 3D double
void xcorr(const double* const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
           const double* const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                 double* const resData,
           eXCorrRes shape = XCORR_RES_FULL);

/***************************************************
 * Normalized cross-correlation
 **************************************************/

// 1D
template <typename T>
void ecorr(T *tplData, size_t tplNrow,
           T *refData, size_t refNrow,
           T *resData,
           eXCorrRes shape = XCORR_RES_FULL);

// 2D
template <typename T>
void ecorr(T *tplData, size_t tplNrow, size_t tplNcol,
           T *refData, size_t refNrow, size_t refNcol,
           T *resData,
           eXCorrRes shape = XCORR_RES_FULL);

// 3D
template <typename T>
void ecorr(T *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
           T *refData, size_t refNrow, size_t refNcol, size_t refNsec,
           T *resData,
           eXCorrRes shape = XCORR_RES_FULL);

// 1D
template <typename T>
void normxcorr(T *tplData, size_t tplNrow,
               T *refData, size_t refNrow,
               T *resData,
               eXCorrRes shape = XCORR_RES_FULL);

// 2D
template <typename T>
void normxcorr(T *tplData, size_t tplNrow, size_t tplNcol,
               T *refData, size_t refNrow, size_t refNcol,
               T *resData,
               eXCorrRes shape = XCORR_RES_FULL);

// 3D
template <typename T>
void normxcorr(T *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
               T *refData, size_t refNrow, size_t refNcol, size_t refNsec,
               T *resData,
               eXCorrRes shape = XCORR_RES_FULL);

/***************************************************
 * Normalized cross-correlation with mask
 **************************************************/

// 1D
template <typename T>
void ecorrm(T *tplData, size_t tplNrow,
            T *refData, size_t refNrow,
            T *resData,
            T *mskData,
            eXCorrRes shape = XCORR_RES_FULL);

// 2D
template <typename T>
void ecorrm(T *tplData, size_t tplNrow, size_t tplNcol,
            T *refData, size_t refNrow, size_t refNcol,
            T *resData,
            T *mskData,
            eXCorrRes shape = XCORR_RES_FULL);

// 3D
template <typename T>
void ecorrm(T *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
            T *refData, size_t refNrow, size_t refNcol, size_t refNsec,
            T *resData,
            T *mskData,
            eXCorrRes shape = XCORR_RES_FULL);

// 1D
template <typename T>
void normxcorrm(T *tplData, size_t tplNrow,
                T *refData, size_t refNrow,
                T *resData,
                T *mskData,
                eXCorrRes shape = XCORR_RES_FULL);

// 2D
template <typename T>
void normxcorrm(T *tplData, size_t tplNrow, size_t tplNcol,
                T *refData, size_t refNrow, size_t refNcol,
                T *resData,
                T *mskData,
                eXCorrRes shape = XCORR_RES_FULL);

// 3D
template <typename T>
void normxcorrm(T *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                T *refData, size_t refNrow, size_t refNcol, size_t refNsec,
                T *resData,
                T *mskData,
                eXCorrRes shape = XCORR_RES_FULL);

/***************************************************
 * Combination for fast computation
 **************************************************/

// 1D float
void normxcorrm_fast(float *tplData, size_t tplNrow,
                     float *refData, size_t refNrow,
                     float *resData,
                     float *mskData);

// 1D double
void normxcorrm_fast(double *tplData, size_t tplNrow,
                     double *refData, size_t refNrow,
                     double *resData,
                     double *mskData);

// 2D float
void normxcorrm_fast(float *tplData, size_t tplNrow, size_t tplNcol,
                     float *refData, size_t refNrow, size_t refNcol,
                     float *resData,
                     float *mskData);

// 2D double
void normxcorrm_fast(double *tplData, size_t tplNrow, size_t tplNcol,
                     double *refData, size_t refNrow, size_t refNcol,
                     double *resData,
                     double *mskData);

// 3D float
void normxcorrm_fast(float *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                     float *refData, size_t refNrow, size_t refNcol, size_t refNsec,
                     float *resData,
                     float *mskData);

// 3D double
void normxcorrm_fast(double *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                     double *refData, size_t refNrow, size_t refNcol, size_t refNsec,
                     double *resData,
                     double *mskData);

/***************************************************
 * Normalized cross-correlation with mask + weight
 **************************************************/

// 1D
template <typename T>
void normxcorrmw(T *tplData, size_t tplNrow, T *tplWght,
                 T *refData, size_t refNrow, T *refWght,
                 T *resData,
                 T *mskData,
                 eXCorrRes shape = XCORR_RES_FULL);

// 2D
template <typename T>
void normxcorrmw(T *tplData, size_t tplNrow, size_t tplNcol, T *tplWght,
                 T *refData, size_t refNrow, size_t refNcol, T *refWght,
                 T *resData,
                 T *mskData,
                 eXCorrRes shape = XCORR_RES_FULL);

// 3D
template <typename T>
void normxcorrmw(T *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec, T *tplWght,
                 T *refData, size_t refNrow, size_t refNcol, size_t refNsec, T *refWght,
                 T *resData,
                 T *mskData,
                 eXCorrRes shape = XCORR_RES_FULL);

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
                         T  overlapRatio);

// 2D
template <typename T>
void regnormxcorrm(const T* const tplData, size_t tplNrow, size_t tplNcol,
                   const T* const refData, size_t refNrow, size_t refNcol,
                   const T* const tplMask,
                   const T* const refMask,
                         T* const resCorr,
                         T* const resOlap,
                         T  overlapRatio);

// 3D
template <typename T>
void regnormxcorrm(const T* const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                   const T* const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                   const T* const tplMask,
                   const T* const refMask,
                         T* const resCorr,
                         T* const resOlap,
                         T  overlapRatio);

} // namespace gem

#endif

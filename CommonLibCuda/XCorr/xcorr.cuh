/***********************************************************************
 *  File:       xcorr.cuh
 *
 *  Purpose:    Header file for cross-correlation functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_DEV_XCORR_CUH__
#define __GEM_DEV_XCORR_CUH__

#include "array.cuh"
#include "xcorr.hpp"

#include <cufft.h>
#include <cuComplex.h>

namespace gem {

void printErrMsgCUFFT(cufftResult_t err);

/***********************************************************************
 * Helper functions for computing correlation
 **********************************************************************/

// modulation
template <typename T1, typename T2>
void cuda_pcorrModulateAndNormalize(const T1* const tplFFT,
                                    const T1* const refFFT,
                                          T1* const resFFT,
                                    size_t fftSize,
                                    T2     nElements);

template <typename T1, typename T2>
void cuda_xcorrModulateAndNormalize(const T1* const tplFFT,
                                    const T1* const refFFT,
                                          T1* const resFFT,
                                    size_t fftSize,
                                    T2     nElements);

// template normalization
template <typename T>
void cuda_xcorrNormTemplateECC(T* const tplData, size_t tplSize);

template <typename T>
void cuda_xcorrNormTemplateECC(      T* const tplData, size_t tplSize,
                               const T* const mskData);

template <typename T>
size_t cuda_xcorrNormTemplateNCC(T* const tplData, size_t tplSize);

template <typename T>
size_t cuda_xcorrNormTemplateNCC(      T* const tplData, size_t tplSize,
                                 const T* const mskData);

template <typename T>
void cuda_xcorrNormTemplateWCC1(T* const tplData, size_t tplSize,
                                T* const tplWght);

template <typename T>
void cuda_xcorrNormTemplateWCC2(T* const tplData, size_t tplSize,
                                T* const tplWght,
                                T* const refWght);

template <typename T>
void cuda_xcorrNormTemplateWCC1(      T* const tplData, size_t tplSize,
                                const T* const mskData,
                                      T* const tplWght);

template <typename T>
void cuda_xcorrNormTemplateWCC2(      T* const tplData, size_t tplSize,
                                const T* const mskData,
                                      T* const tplWght,
                                      T* const refWght);

// combine correlations
template <typename T>
void cuda_xcorrCombineResultCCC(      T* const resData,
                                const T* const corrData1, T weight1,
                                const T* const corrData2, T weight2,
                                const T* const corrData3, T weight3,
                                size_t resSize);

template <typename T>
void cuda_xcorrCombineResultECC(      T* const resData,
                                const T* const fftTmp1,
                                size_t resSize);

template <typename T>
void cuda_xcorrCombineResultNCC(      T* const resData,
                                const T* const fftTmp1,
                                const T* const fftTmp2,
                                size_t mskArea,
                                size_t resSize);

template <typename T>
void cuda_xcorrCombineResultWCC(      T* const resData,
                                const T* const fftTmp1,
                                const T* const fftTmp2,
                                T      tplSum,
                                size_t resSize);

template <typename T>
void cuda_xcorrCombineResultNReg(      T* const resCorr,
                                 const T* const resOlap,
                                 const T* const numTpl,
                                 const T* const numRef,
                                 const T* const denTpl,
                                 const T* const denRef,
                                 size_t resSize,
                                 T      overlapRatio);

// merge result
template <typename T>
void cuda_xcorrMergeResult(const T*      const resData,
                                 T*      const resDataAbsMax,
                                 size_t* const resDataMaxInd,
                           size_t        resSize,
                           size_t        indx,
                           eXCorrMerge   bAbs);

/***********************************************************************
 * Phase-correlation
 **********************************************************************/

// size(tplData) = size(refData) = size(resData) = refNrow
void cuda_pcorr(const float*  const tplData,
                const float*  const refData,
                      float*  const resData,
                size_t refNrow);

void cuda_pcorr(const double* const tplData,
                const double* const refData,
                      double* const resData,
                size_t refNrow);

// size(tplData) = size(refData) = size(resData) = [refNrow refNcol]
void cuda_pcorr(const float*  const tplData,
                const float*  const refData,
                      float*  const resData,
                size_t refNrow, size_t refNcol);

void cuda_pcorr(const double* const tplData,
                const double* const refData,
                      double* const resData,
                size_t refNrow, size_t refNcol);

// size(tplData) = size(refData) = size(resData) = [refNrow refNcol refNsec]
void cuda_pcorr(const float*  const tplData,
                const float*  const refData,
                      float*  const resData,
                size_t refNrow, size_t refNcol, size_t refNsec);

void cuda_pcorr(const double* const tplData,
                const double* const refData,
                      double* const resData,
                size_t refNrow, size_t refNcol, size_t refNsec);

/***********************************************************************
 * Direct correlation
 **********************************************************************/

// xcorr
template <typename T>
T cuda_xcorr_direct(const T* const tplData,
                    const T* const refData,
                    size_t length);

// ecorr
template <typename T>
T cuda_ecorr_direct(const T* const tplData,
                    const T* const refData,
                    size_t length);

// ecorrm
template <typename T>
T cuda_ecorrm_direct(const T* const tplData,
                     const T* const refData,
                     const T* const mskData,
                     size_t length);

// normxcorr
template <typename T>
T cuda_normxcorr_direct(const T* const tplData,
                        const T* const refData,
                        size_t length);

// normxcorrm
template <typename T>
T cuda_normxcorrm_direct(const T* const tplData,
                         const T* const refData,
                         const T* const mskData,
                         size_t length);

// normxcorrmw
template <typename T>
T cuda_normxcorrmw_direct(const T* const tplData, const T* const tplWght,
                          const T* const refData, const T* const refWght,
                          const T* const mskData,
                          size_t length);

/***********************************************************************
 * Cross-correlation
 **********************************************************************/

// 1D float
void cuda_xcorr(const float*  const tplData, size_t tplNrow,
                const float*  const refData, size_t refNrow,
                      float*  const resData,
                eXCorrRes shape = XCORR_RES_FULL);

// 1D double
void cuda_xcorr(const double* const tplData, size_t tplNrow,
                const double* const refData, size_t refNrow,
                      double* const resData,
                eXCorrRes shape = XCORR_RES_FULL);

// 2D float
void cuda_xcorr(const float*  const tplData, size_t tplNrow, size_t tplNcol,
                const float*  const refData, size_t refNrow, size_t refNcol,
                      float*  const resData,
                eXCorrRes shape = XCORR_RES_FULL);

// 2D double
void cuda_xcorr(const double* const tplData, size_t tplNrow, size_t tplNcol,
                const double* const refData, size_t refNrow, size_t refNcol,
                      double* const resData,
                eXCorrRes shape = XCORR_RES_FULL);

// 3D float
void cuda_xcorr(const float*  const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                const float*  const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                      float*  const resData,
                eXCorrRes shape = XCORR_RES_FULL);

// 3D double
void cuda_xcorr(const double* const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                const double* const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                      double* const resData,
                eXCorrRes shape = XCORR_RES_FULL);

/***********************************************************************
 * Normalized cross-correlation
 **********************************************************************/

// 1D
template <typename T>
void cuda_ecorr(T *dev_tplData, size_t tplNrow,
                T *dev_refData, size_t refNrow,
                T *dev_resData,
                eXCorrRes shape = XCORR_RES_FULL);

// 2D
template <typename T>
void cuda_ecorr(T *dev_tplData, size_t tplNrow, size_t tplNcol,
                T *dev_refData, size_t refNrow, size_t refNcol,
                T *dev_resData,
                eXCorrRes shape = XCORR_RES_FULL);

// 3D
template <typename T>
void cuda_ecorr(T *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                T *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                T *dev_resData,
                eXCorrRes shape = XCORR_RES_FULL);

// 1D
template <typename T>
void cuda_normxcorr(T *dev_tplData, size_t tplNrow,
                    T *dev_refData, size_t refNrow,
                    T *dev_resData,
                    eXCorrRes shape = XCORR_RES_FULL);

// 2D
template <typename T>
void cuda_normxcorr(T *dev_tplData, size_t tplNrow, size_t tplNcol,
                    T *dev_refData, size_t refNrow, size_t refNcol,
                    T *dev_resData,
                    eXCorrRes shape = XCORR_RES_FULL);

// 3D
template <typename T>
void cuda_normxcorr(T *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                    T *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                    T *dev_resData,
                    eXCorrRes shape = XCORR_RES_FULL);

/***********************************************************************
 * Normalized cross-correlation with mask
 **********************************************************************/

// 1D
template <typename T>
void cuda_ecorrm(T *dev_tplData, size_t tplNrow,
                 T *dev_refData, size_t refNrow,
                 T *dev_resData,
                 T *dev_mskData,
                 eXCorrRes shape = XCORR_RES_FULL);

// 2D
template <typename T>
void cuda_ecorrm(T *dev_tplData, size_t tplNrow, size_t tplNcol,
                 T *dev_refData, size_t refNrow, size_t refNcol,
                 T *dev_resData,
                 T *dev_mskData,
                 eXCorrRes shape = XCORR_RES_FULL);

// 3D
template <typename T>
void cuda_ecorrm(T *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                 T *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                 T *dev_resData,
                 T *dev_mskData,
                 eXCorrRes shape = XCORR_RES_FULL);

// 1D
template <typename T>
void cuda_normxcorrm(T *dev_tplData, size_t tplNrow,
                     T *dev_refData, size_t refNrow,
                     T *dev_resData,
                     T *dev_mskData,
                     eXCorrRes shape = XCORR_RES_FULL);

// 2D
template <typename T>
void cuda_normxcorrm(T *dev_tplData, size_t tplNrow, size_t tplNcol,
                     T *dev_refData, size_t refNrow, size_t refNcol,
                     T *dev_resData,
                     T *dev_mskData,
                     eXCorrRes shape = XCORR_RES_FULL);

// 3D
template <typename T>
void cuda_normxcorrm(T *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                     T *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                     T *dev_resData,
                     T *dev_mskData,
                     eXCorrRes shape = XCORR_RES_FULL);

/***********************************************************************
 * Combination for fast computation
 **********************************************************************/

// 1D float
void cuda_normxcorrm_fast(float *dev_tplData, size_t tplNrow,
                          float *dev_refData, size_t refNrow,
                          float *dev_resData,
                          float *dev_mskData);

// 1D double
void cuda_normxcorrm_fast(double *dev_tplData, size_t tplNrow,
                          double *dev_refData, size_t refNrow,
                          double *dev_resData,
                          double *dev_mskData);

// 2D float
void cuda_normxcorrm_fast(float *dev_tplData, size_t tplNrow, size_t tplNcol,
                          float *dev_refData, size_t refNrow, size_t refNcol,
                          float *dev_resData,
                          float *dev_mskData);

// 2D double
void cuda_normxcorrm_fast(double *dev_tplData, size_t tplNrow, size_t tplNcol,
                          double *dev_refData, size_t refNrow, size_t refNcol,
                          double *dev_resData,
                          double *dev_mskData);

// 3D float
void cuda_normxcorrm_fast(float *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                          float *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                          float *dev_resData,
                          float *dev_mskData);

// 3D double
void cuda_normxcorrm_fast(double *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                          double *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                          double *dev_resData,
                          double *dev_mskData);

/***********************************************************************
 * Normalized cross-correlation with mask + weight
 **********************************************************************/

// 1D
template <typename T>
void cuda_normxcorrmw(T *dev_tplData, size_t tplNrow, T *dev_tplWght,
                      T *dev_refData, size_t refNrow, T *dev_refWght,
                      T *dev_resData,
                      T *dev_mskData,
                      eXCorrRes shape = XCORR_RES_FULL);

// 2D
template <typename T>
void cuda_normxcorrmw(T *dev_tplData, size_t tplNrow, size_t tplNcol, T *dev_tplWght,
                      T *dev_refData, size_t refNrow, size_t refNcol, T *dev_refWght,
                      T *dev_resData,
                      T *dev_mskData,
                      eXCorrRes shape = XCORR_RES_FULL);

// 3D
template <typename T>
void cuda_normxcorrmw(T *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec, T *dev_tplWght,
                      T *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec, T *dev_refWght,
                      T *dev_resData,
                      T *dev_mskData,
                      eXCorrRes shape = XCORR_RES_FULL);

/***************************************************
 * Masked registration using NCC
 **************************************************/

// 1D
template <typename T>
void cuda_regnormxcorrm(const T* const tplData, size_t tplNrow,
                        const T* const refData, size_t refNrow,
                        const T* const tplMask,
                        const T* const refMask,
                              T* const resCorr,
                              T* const resOlap,
                              T  overlapRatio);

// 2D
template <typename T>
void cuda_regnormxcorrm(const T* const tplData, size_t tplNrow, size_t tplNcol,
                        const T* const refData, size_t refNrow, size_t refNcol,
                        const T* const tplMask,
                        const T* const refMask,
                              T* const resCorr,
                              T* const resOlap,
                              T  overlapRatio);

// 3D
template <typename T>
void cuda_regnormxcorrm(const T* const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                        const T* const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                        const T* const tplMask,
                        const T* const refMask,
                              T* const resCorr,
                              T* const resOlap,
                              T  overlapRatio);

} // namespace gem

#endif

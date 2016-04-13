/***********************************************************************
 *  File:       xcorr_pcorr.cu
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

} // namespace gem

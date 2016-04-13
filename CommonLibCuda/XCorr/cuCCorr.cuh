/***********************************************************************
 *  File:       cuCCorr.cuh
 *
 *  Purpose:    Header file for CCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUCCORR_CUH__
#define __GEM_CUCCORR_CUH__

#include "xcorr.cuh"
#include "cuWCorr.cuh"
#include "cuXCorr.cuh"

namespace gem {

class cuCCorrSingle
{
private:
    size_t          refSize;

    float           *dev_resDataPadded;

    float           *dev_resDataAbsMaxPadded;
    size_t          *dev_resDataMaxIndPadded;

    float           *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    cuWCorrSingle   *corrDensity;
    cuXCorrSingle   *corrSurface;
    cuXCorrSingle   *corrPenalty;

    float           wDensity;
    float           wSurface;
    float           wPenalty;

public:
     cuCCorrSingle(size_t refNrowInput,
                   size_t refNcolInput,
                   size_t refNsecInput);
    ~cuCCorrSingle() {};

    void memAlloc(void);
    void memFree (void);

    void setSizeTpl(size_t tplNrowInput,
                    size_t tplNcolInput,
                    size_t tplNsecInput);

    void rotateTpl(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
    void normalizeTpl(void);
    void prepareTpl(void);
    void computeCorr(void);

    void mergeResult(size_t indx, eXCorrMerge bAbs);
    void mergeResultGlobal(float  *resDataAbsMaxPaddedGlobal,
                           size_t *resDataMaxIndPaddedGlobal,
                           eXCorrMerge bAbs);

    void setRotTexture(bool bRotTextureInput);
    void setWeight(float wDensityInput, float wSurfaceInput, float wPenaltyInput);

    cuWCorrSingle*  getRefCorrDensity(void) { return corrDensity; };
    cuXCorrSingle*  getRefCorrSurface(void) { return corrSurface; };
    cuXCorrSingle*  getRefCorrPenalty(void) { return corrPenalty; };
};

} // namespace gem

#endif

/***********************************************************************
 *  File:       cCCorr.hpp
 *
 *  Purpose:    Header file for CCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CCCORR_HPP__
#define __GEM_CCCORR_HPP__

#include "xcorr.hpp"
#include "cWCorr.hpp"
#include "cXCorr.hpp"

namespace gem {

class cCCorrSingle
{
private:
    size_t          refSize;

    float           *resDataPadded;

    float           *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    cWCorrSingle   *corrDensity;
    cXCorrSingle   *corrSurface;
    cXCorrSingle   *corrPenalty;

    float           wDensity;
    float           wSurface;
    float           wPenalty;

public:
     cCCorrSingle(size_t refNrowInput,
                  size_t refNcolInput,
                  size_t refNsecInput);
    ~cCCorrSingle() {};

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

    void setWeight(float wDensityInput, float wSurfaceInput, float wPenaltyInput);

    cWCorrSingle*  getRefCorrDensity(void) { return corrDensity; };
    cXCorrSingle*  getRefCorrSurface(void) { return corrSurface; };
    cXCorrSingle*  getRefCorrPenalty(void) { return corrPenalty; };
};

} // namespace gem

#endif

/***********************************************************************
 *  File:       cCCorr.cpp
 *
 *  Purpose:    Implementation of CCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cCCorr.hpp"

namespace gem {

cCCorrSingle::cCCorrSingle(size_t refNrowInput,
                            size_t refNcolInput,
                            size_t refNsecInput)
{
    refSize = refNrowInput * refNcolInput * refNsecInput;

    resDataPadded = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    corrDensity = new cWCorrSingle(refNrowInput, refNcolInput, refNsecInput);
    corrSurface = new cXCorrSingle(refNrowInput, refNcolInput, refNsecInput);
    corrPenalty = new cXCorrSingle(refNrowInput, refNcolInput, refNsecInput);

    wDensity = 0.5f;
    wSurface = 0.5f;
    wPenalty = 0.5f;
}

void cCCorrSingle::memAlloc(void)
{
    array_new(resDataPadded, refSize);

    array_new_zero(resDataAbsMaxPadded, refSize);
    array_new     (resDataMaxIndPadded, refSize);

    corrDensity->memAlloc();
    corrSurface->memAlloc();
    corrPenalty->memAlloc();
}

void cCCorrSingle::memFree(void)
{
    array_delete(resDataPadded);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    corrDensity->memFree();
    corrSurface->memFree();
    corrPenalty->memFree();

    delete corrDensity;
    delete corrSurface;
    delete corrPenalty;
}

void cCCorrSingle::setSizeTpl(size_t tplNrowInput,
                              size_t tplNcolInput,
                              size_t tplNsecInput)
{
    corrDensity->setSizeTpl(tplNrowInput, tplNcolInput, tplNsecInput);
    corrSurface->setSizeTpl(tplNrowInput, tplNcolInput, tplNsecInput);
    corrPenalty->setSizeTpl(tplNrowInput, tplNcolInput, tplNsecInput);
}

void cCCorrSingle::rotateTpl(float alpha, float beta, float gamma, eInter inter)
{
    corrDensity->tplCreateByRot3D(alpha, beta, gamma, inter);
    corrDensity->mskCreateByRot3D(alpha, beta, gamma, inter);
    corrDensity->wghtTplCreateByRot3D(alpha, beta, gamma, inter);
    corrDensity->wghtRefCreateByRot3D(alpha, beta, gamma, inter);

    corrSurface->tplCreateByRot3D(alpha, beta, gamma, inter);
    corrPenalty->tplCreateByRot3D(alpha, beta, gamma, inter);
}

void cCCorrSingle::normalizeTpl(void)
{
    corrDensity->normalizeTpl();
}

void cCCorrSingle::prepareTpl(void)
{
    corrDensity->prepareTplRefWght();
    corrSurface->prepareTpl();
    corrPenalty->prepareTpl();
}

void cCCorrSingle::computeCorr(void)
{
    corrDensity->computeCorr();
    corrSurface->computeCorr();
    corrPenalty->computeCorr();
}

void cCCorrSingle::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    xcorrCombineResultCCC(resDataPadded,
                          corrDensity->getResultRef(), wDensity,
                          corrSurface->getResultRef(), wSurface,
                          corrPenalty->getResultRef(), wPenalty,
                          refSize);

    xcorrMergeResult(resDataPadded,
                     resDataAbsMaxPadded,
                     resDataMaxIndPadded,
                     refSize,
                     indx,
                     bAbs);
}

void cCCorrSingle::mergeResultGlobal(float  *resDataAbsMaxPaddedGlobal,
                                     size_t *resDataMaxIndPaddedGlobal,
                                     eXCorrMerge bAbs)
{
    xcorrMergeResultGlobal(resDataAbsMaxPadded,
                           resDataMaxIndPadded,
                           resDataAbsMaxPaddedGlobal,
                           resDataMaxIndPaddedGlobal,
                           refSize,
                           bAbs);
}

void cCCorrSingle::setWeight(float wDensityInput, float wSurfaceInput, float wPenaltyInput)
{
    wDensity = wDensityInput;
    wSurface = wSurfaceInput;
    wPenalty = wPenaltyInput;
}

} // namespace gem

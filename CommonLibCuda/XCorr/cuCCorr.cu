/***********************************************************************
 *  File:       cuCCorr.cu
 *
 *  Purpose:    Implementation of CCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuCCorr.cuh"

namespace gem {

cuCCorrSingle::cuCCorrSingle(size_t refNrowInput,
                             size_t refNcolInput,
                             size_t refNsecInput)
{
    refSize = refNrowInput * refNcolInput * refNsecInput;

    dev_resDataPadded = NULL;

    dev_resDataAbsMaxPadded = NULL;
    dev_resDataMaxIndPadded = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    corrDensity = new cuWCorrSingle(refNrowInput, refNcolInput, refNsecInput);
    corrSurface = new cuXCorrSingle(refNrowInput, refNcolInput, refNsecInput);
    corrPenalty = new cuXCorrSingle(refNrowInput, refNcolInput, refNsecInput);

    wDensity = 0.5f;
    wSurface = 0.5f;
    wPenalty = 0.5f;
}

void cuCCorrSingle::memAlloc(void)
{
    cuda_arrayDev_new(dev_resDataPadded, refSize);

    cuda_arrayDev_new_zero(dev_resDataAbsMaxPadded, refSize);
    cuda_arrayDev_new     (dev_resDataMaxIndPadded, refSize);

    array_new(resDataAbsMaxPadded, refSize);
    array_new(resDataMaxIndPadded, refSize);

    corrDensity->memAlloc();
    corrSurface->memAlloc();
    corrPenalty->memAlloc();
}

void cuCCorrSingle::memFree(void)
{
    cuda_arrayDev_delete(dev_resDataPadded);

    cuda_arrayDev_delete(dev_resDataAbsMaxPadded);
    cuda_arrayDev_delete(dev_resDataMaxIndPadded);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    corrDensity->memFree();
    corrSurface->memFree();
    corrPenalty->memFree();

    delete corrDensity;
    delete corrSurface;
    delete corrPenalty;
}

void cuCCorrSingle::setSizeTpl(size_t tplNrowInput,
                               size_t tplNcolInput,
                               size_t tplNsecInput)
{
    corrDensity->setSizeTpl(tplNrowInput, tplNcolInput, tplNsecInput);
    corrSurface->setSizeTpl(tplNrowInput, tplNcolInput, tplNsecInput);
    corrPenalty->setSizeTpl(tplNrowInput, tplNcolInput, tplNsecInput);
}

void cuCCorrSingle::rotateTpl(float alpha, float beta, float gamma, eInter inter)
{
    corrDensity->tplCreateByRot3D(alpha, beta, gamma, inter);
    corrDensity->mskCreateByRot3D(alpha, beta, gamma, inter);
    corrDensity->wghtTplCreateByRot3D(alpha, beta, gamma, inter);
    corrDensity->wghtRefCreateByRot3D(alpha, beta, gamma, inter);

    corrSurface->tplCreateByRot3D(alpha, beta, gamma, inter);
    corrPenalty->tplCreateByRot3D(alpha, beta, gamma, inter);
}

void cuCCorrSingle::normalizeTpl(void)
{
    corrDensity->normalizeTpl();
}

void cuCCorrSingle::prepareTpl(void)
{
    corrDensity->prepareTplRefWght();
    corrSurface->prepareTpl();
    corrPenalty->prepareTpl();
}

void cuCCorrSingle::computeCorr(void)
{
    corrDensity->computeCorr();
    corrSurface->computeCorr();
    corrPenalty->computeCorr();
}

void cuCCorrSingle::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    cuda_xcorrCombineResultCCC(dev_resDataPadded,
                               corrDensity->getResultRef(), wDensity,
                               corrSurface->getResultRef(), wSurface,
                               corrPenalty->getResultRef(), wPenalty,
                               refSize);

    cuda_xcorrMergeResult(dev_resDataPadded,
                          dev_resDataAbsMaxPadded,
                          dev_resDataMaxIndPadded,
                          refSize,
                          indx,
                          bAbs);
}

void cuCCorrSingle::mergeResultGlobal(float  *resDataAbsMaxPaddedGlobal,
                                      size_t *resDataMaxIndPaddedGlobal,
                                      eXCorrMerge bAbs)
{
    cuda_array_memcpy_d2h(resDataAbsMaxPadded, dev_resDataAbsMaxPadded, refSize);
    cuda_array_memcpy_d2h(resDataMaxIndPadded, dev_resDataMaxIndPadded, refSize);

    xcorrMergeResultGlobal(resDataAbsMaxPadded,
                           resDataMaxIndPadded,
                           resDataAbsMaxPaddedGlobal,
                           resDataMaxIndPaddedGlobal,
                           refSize,
                           bAbs);
}

void cuCCorrSingle::setRotTexture(bool bRotTextureInput)
{
    corrDensity->setRotTexture(bRotTextureInput);
    corrSurface->setRotTexture(bRotTextureInput);
    corrPenalty->setRotTexture(bRotTextureInput);
}

void cuCCorrSingle::setWeight(float wDensityInput, float wSurfaceInput, float wPenaltyInput)
{
    wDensity = wDensityInput;
    wSurface = wSurfaceInput;
    wPenalty = wPenaltyInput;
}

} // namespace gem

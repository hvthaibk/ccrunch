/***********************************************************************
 *  File:       cRegCorr.hpp
 *
 *  Purpose:    Header file for a correlation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CREG_CORR_HPP__
#define __GEM_CREG_CORR_HPP__

#include "xcorr.hpp"
#include "cData3XCorr.hpp"

#ifdef __GEM_USE_CUDA__
#include "xcorr.cuh"
#include "cuData3.hpp"
#endif

namespace gem {

template <typename T>
class cRegCorr
{
private:

public:
     cRegCorr(void) {};
    ~cRegCorr()     {};

    T       computeXCC(const cData3XCorr<T>& tpl, const cData3XCorr<T>& ref);
    void    computeXCC(const cData3XCorr<T>& tpl, const cData3XCorr<T>& ref, cData3XCorr<T>& res, eXCorrRes shape = XCORR_RES_FULL);

    T       computeECC(const cData3XCorr<T>& tpl, const cData3XCorr<T>& ref,                                                        bool withMask = true);
    void    computeECC(const cData3XCorr<T>& tpl, const cData3XCorr<T>& ref, cData3XCorr<T>& res, eXCorrRes shape = XCORR_RES_FULL, bool withMask = true);

    T       computeNCC(const cData3XCorr<T>& tpl, const cData3XCorr<T>& ref,                                                        bool withMask = true);
    void    computeNCC(const cData3XCorr<T>& tpl, const cData3XCorr<T>& ref, cData3XCorr<T>& res, eXCorrRes shape = XCORR_RES_FULL, bool withMask = true);

    T       computeWCC(const cData3XCorr<T>& tpl, const cData3XCorr<T>& ref);
    void    computeWCC(const cData3XCorr<T>& tpl, const cData3XCorr<T>& ref, cData3XCorr<T>& res, eXCorrRes shape = XCORR_RES_FULL);

    void    computePCC(const cData3XCorr<T>& tpl, const cData3XCorr<T>& ref, cData3XCorr<T>& res);

    void    computeREG(const cData3XCorr<T>& tpl, const cData3XCorr<T>& ref, cData3XCorr<T>& res, T overlapRatio);

#ifdef __GEM_USE_CUDA__

    T       computeXCC(const cuData3<T>& tpl, const cuData3<T>& ref);
    void    computeXCC(const cuData3<T>& tpl, const cuData3<T>& ref, cuData3<T>& res, eXCorrRes shape = XCORR_RES_FULL);

    T       computeECC(const cuData3<T>& tpl, const cuData3<T>& ref);
    T       computeECC(const cuData3<T>& tpl, const cuData3<T>& ref, const cuData3<T>& msk);
    void    computeECC(const cuData3<T>& tpl, const cuData3<T>& ref,
                             cuData3<T>& res, eXCorrRes shape = XCORR_RES_FULL);
    void    computeECC(const cuData3<T>& tpl, const cuData3<T>& ref, const cuData3<T>& msk,
                             cuData3<T>& res, eXCorrRes shape = XCORR_RES_FULL);

    T       computeNCC(const cuData3<T>& tpl, const cuData3<T>& ref);
    T       computeNCC(const cuData3<T>& tpl, const cuData3<T>& ref, const cuData3<T>& msk);
    void    computeNCC(const cuData3<T>& tpl, const cuData3<T>& ref,
                             cuData3<T>& res, eXCorrRes shape = XCORR_RES_FULL);
    void    computeNCC(const cuData3<T>& tpl, const cuData3<T>& ref, const cuData3<T>& msk,
                             cuData3<T>& res, eXCorrRes shape = XCORR_RES_FULL);

    T       computeWCC(const cuData3<T>& tpl, const cuData3<T>& ref, const cuData3<T>& msk, const cuData3<T>& tplWght, const cuData3<T>& refWght);
    void    computeWCC(const cuData3<T>& tpl, const cuData3<T>& ref, const cuData3<T>& msk, const cuData3<T>& tplWght, const cuData3<T>& refWght,
                             cuData3<T>& res, eXCorrRes shape = XCORR_RES_FULL);

    void    computePCC(const cuData3<T>& tpl, const cuData3<T>& ref, cuData3<T>& res);

    void    computeREG(const cuData3<T>& tpl, const cuData3<T>& ref, cuData3<T>& res, T overlapRatio);

#endif
};

} // namespace gem

#endif

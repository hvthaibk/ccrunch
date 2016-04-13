/***********************************************************************
 *  File:       cRegCorr.cpp
 *
 *  Purpose:    Implementation of a correlation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cRegCorr.hpp"

namespace gem {

template <typename T>
T cRegCorr<T>::computeXCC(const cData3XCorr<T>& tpl,
                          const cData3XCorr<T>& ref)
{
    const std::string   funcName("T cRegCorr<T>::computeXCC("
                                    "const cData3XCorr<T>& tpl, "
                                    "const cData3XCorr<T>& ref)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    return xcorr_direct(tpl.getAddrData(),
                        ref.getAddrData(),
                        tpl.getNelement());
}

template <typename T>
void cRegCorr<T>::computeXCC(const cData3XCorr<T>& tpl,
                             const cData3XCorr<T>& ref,
                                   cData3XCorr<T>& res,
                             eXCorrRes shape)
{
    const std::string   funcName("void cRegCorr<T>::computeXCC("
                                    "const cData3XCorr<T>& tpl, "
                                    "const cData3XCorr<T>& ref, "
                                    "cData3XCorr<T>& res, "
                                    "eXCorrRes shape)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() <= ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    switch (shape) {
        case XCORR_RES_FULL:
            res.memReAlloc(ref.getSize() + tpl.getSize() - cVector3<size_t>(1,1,1));
            break;
        case XCORR_RES_VALID:
            res.memReAlloc(ref.getSize() - tpl.getSize() + cVector3<size_t>(1,1,1));
            break;
        default:
            ERROR(funcName, "unsupported xcorr mode");
    }

    switch (tpl.getDimension()) {
        case 1:
            xcorr(tpl.getAddrData(), tpl.getNrow(),
                  ref.getAddrData(), ref.getNrow(),
                  res.getAddrData(),
                  shape);
            break;
        case 2:
            xcorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(),
                  ref.getAddrData(), ref.getNrow(), ref.getNcol(),
                  res.getAddrData(),
                  shape);
            break;
        case 3:
            xcorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                  ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                  res.getAddrData(),
                  shape);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
T cRegCorr<T>::computeECC(const cData3XCorr<T>& tpl,
                          const cData3XCorr<T>& ref,
                          bool withMask)
{
    const std::string   funcName("T cRegCorr<T>::computeECC("
                                    "const cData3XCorr<T>& tpl, "
                                    "const cData3XCorr<T>& ref, "
                                    "bool withMask)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    if (withMask) {
        require(tpl.getSize() == tpl.maskGetSize(),
                funcName + ": tpl and its mask have different sizes");

        return ecorrm_direct(tpl.getAddrData(),
                             ref.getAddrData(),
                             tpl.maskGetAddrData(),
                             tpl.getNelement());
    }
    else {
        return ecorr_direct(tpl.getAddrData(),
                            ref.getAddrData(),
                            tpl.getNelement());
    }
}

template <typename T>
void cRegCorr<T>::computeECC(const cData3XCorr<T>& tpl,
                             const cData3XCorr<T>& ref,
                                   cData3XCorr<T>& res,
                             eXCorrRes shape,
                             bool withMask)
{
    const std::string   funcName("void cRegCorr<T>::computeECC("
                                    "const cData3XCorr<T>& tpl, "
                                    "const cData3XCorr<T>& ref, "
                                    "cData3XCorr<T>& res, "
                                    "eXCorrRes shape, "
                                    "bool withMask)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() <= ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    switch (shape) {
        case XCORR_RES_FULL:
            res.memReAlloc(ref.getSize() + tpl.getSize() - cVector3<size_t>(1,1,1));
            break;
        case XCORR_RES_VALID:
            res.memReAlloc(ref.getSize() - tpl.getSize() + cVector3<size_t>(1,1,1));
            break;
        default:
            ERROR(funcName, "unsupported xcorr mode");
    }

    if (withMask) {
        require(tpl.getSize() == tpl.maskGetSize(),
                funcName + ": tpl and its mask have different sizes");

        switch (tpl.getDimension()) {
            case 1:
                ecorrm(tpl.getAddrData(), tpl.getNrow(),
                       ref.getAddrData(), ref.getNrow(),
                       res.getAddrData(),
                       tpl.maskGetAddrData(),
                       shape);
                break;
            case 2:
                ecorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(),
                       ref.getAddrData(), ref.getNrow(), ref.getNcol(),
                       res.getAddrData(),
                       tpl.maskGetAddrData(),
                       shape);
                break;
            case 3:
                ecorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                       ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                       res.getAddrData(),
                       tpl.maskGetAddrData(),
                       shape);
                break;
            default:
                ERROR(funcName, "unsupported dimension");
        }
    }
    else {
        switch (tpl.getDimension()) {
            case 1:
                ecorr(tpl.getAddrData(), tpl.getNrow(),
                      ref.getAddrData(), ref.getNrow(),
                      res.getAddrData(),
                      shape);
                break;
            case 2:
                ecorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(),
                      ref.getAddrData(), ref.getNrow(), ref.getNcol(),
                      res.getAddrData(),
                      shape);
                break;
            case 3:
                ecorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                      ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                      res.getAddrData(),
                      shape);
                break;
            default:
                ERROR(funcName, "unsupported dimension");
        }
    }
}

template <typename T>
T cRegCorr<T>::computeNCC(const cData3XCorr<T>& tpl,
                          const cData3XCorr<T>& ref,
                          bool withMask)
{
    const std::string   funcName("T cRegCorr<T>::computeNCC("
                                    "const cData3XCorr<T>& tpl, "
                                    "const cData3XCorr<T>& ref, "
                                    "bool withMask)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    if (withMask) {
        require(tpl.getSize() == tpl.maskGetSize(),
                funcName + ": tpl and its mask have different sizes");

        return normxcorrm_direct(tpl.getAddrData(),
                                 ref.getAddrData(),
                                 tpl.maskGetAddrData(),
                                 tpl.getNelement());
    }
    else {
        return normxcorr_direct(tpl.getAddrData(),
                                ref.getAddrData(),
                                tpl.getNelement());
    }
}

template <typename T>
void cRegCorr<T>::computeNCC(const cData3XCorr<T>& tpl,
                             const cData3XCorr<T>& ref,
                                   cData3XCorr<T>& res,
                             eXCorrRes shape,
                             bool withMask)
{
    const std::string   funcName("void cRegCorr<T>::computeNCC("
                                    "const cData3XCorr<T>& tpl, "
                                    "const cData3XCorr<T>& ref, "
                                    "cData3XCorr<T>& res, "
                                    "eXCorrRes shape, "
                                    "bool withMask)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() <= ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    switch (shape) {
        case XCORR_RES_FULL:
            res.memReAlloc(ref.getSize() + tpl.getSize() - cVector3<size_t>(1,1,1));
            break;
        case XCORR_RES_VALID:
            res.memReAlloc(ref.getSize() - tpl.getSize() + cVector3<size_t>(1,1,1));
            break;
        default:
            ERROR(funcName, "unsupported xcorr mode");
    }

    if (withMask) {
        require(tpl.getSize() == tpl.maskGetSize(),
                funcName + ": tpl and its mask have different sizes");

        switch (tpl.getDimension()) {
            case 1:
                normxcorrm(tpl.getAddrData(), tpl.getNrow(),
                           ref.getAddrData(), ref.getNrow(),
                           res.getAddrData(),
                           tpl.maskGetAddrData(),
                           shape);
                break;
            case 2:
                normxcorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(),
                           ref.getAddrData(), ref.getNrow(), ref.getNcol(),
                           res.getAddrData(),
                           tpl.maskGetAddrData(),
                           shape);
                break;
            case 3:
                normxcorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                           ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                           res.getAddrData(),
                           tpl.maskGetAddrData(),
                           shape);
                break;
            default:
                ERROR(funcName, "unsupported dimension");
        }
    }
    else {
        switch (tpl.getDimension()) {
            case 1:
                normxcorr(tpl.getAddrData(), tpl.getNrow(),
                          ref.getAddrData(), ref.getNrow(),
                          res.getAddrData(),
                          shape);
                break;
            case 2:
                normxcorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(),
                          ref.getAddrData(), ref.getNrow(), ref.getNcol(),
                          res.getAddrData(),
                          shape);
                break;
            case 3:
                normxcorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                          ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                          res.getAddrData(),
                          shape);
                break;
            default:
                ERROR(funcName, "unsupported dimension");
        }
    }
}

template <typename T>
T cRegCorr<T>::computeWCC(const cData3XCorr<T>& tpl,
                          const cData3XCorr<T>& ref)
{
    const std::string   funcName("T cRegCorr<T>::computeWCC("
                                    "const cData3XCorr<T>& tpl, "
                                    "const cData3XCorr<T>& ref)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");
    require(tpl.getSize() == tpl.maskGetSize(),
            funcName + ": tpl and its mask have different sizes");
    require(tpl.getSize() == tpl.wghtTplGetSize(),
            funcName + ": tpl and its wghtTpl have different sizes");
    require(tpl.getSize() == tpl.wghtRefGetSize(),
            funcName + ": tpl and its wghtRef have different sizes");

    return normxcorrmw_direct(tpl.getAddrData(), tpl.wghtTplGetAddrData(),
                              ref.getAddrData(), tpl.wghtRefGetAddrData(),
                              tpl.maskGetAddrData(),
                              tpl.getNelement());
}

template <typename T>
void cRegCorr<T>::computeWCC(const cData3XCorr<T>& tpl,
                             const cData3XCorr<T>& ref,
                                   cData3XCorr<T>& res,
                             eXCorrRes shape)
{
    const std::string   funcName("void cRegCorr<T>::computeWCC("
                                    "const cData3XCorr<T>& tpl, "
                                    "const cData3XCorr<T>& ref, "
                                    "cData3XCorr<T>& res, "
                                    "eXCorrRes shape)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() <= ref.getSize(),
            funcName + ": tpl and ref have different sizes");
    require(tpl.getSize() == tpl.maskGetSize(),
            funcName + ": tpl and its mask have different sizes");
    require(tpl.getSize() == tpl.wghtTplGetSize(),
            funcName + ": tpl and its wghtTpl have different sizes");
    require(tpl.getSize() == tpl.wghtRefGetSize(),
            funcName + ": tpl and its wghtRef have different sizes");

    switch (shape) {
        case XCORR_RES_FULL:
            res.memReAlloc(ref.getSize() + tpl.getSize() - cVector3<size_t>(1,1,1));
            break;
        case XCORR_RES_VALID:
            res.memReAlloc(ref.getSize() - tpl.getSize() + cVector3<size_t>(1,1,1));
            break;
        default:
            ERROR(funcName, "unsupported xcorr mode");
    }

    switch (tpl.getDimension()) {
        case 1:
            normxcorrmw(tpl.getAddrData(), tpl.getNrow(), tpl.wghtTplGetAddrData(),
                        ref.getAddrData(), ref.getNrow(), tpl.wghtRefGetAddrData(),
                        res.getAddrData(),
                        tpl.maskGetAddrData(),
                        shape);
            break;
        case 2:
            normxcorrmw(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.wghtTplGetAddrData(),
                        ref.getAddrData(), ref.getNrow(), ref.getNcol(), tpl.wghtRefGetAddrData(),
                        res.getAddrData(),
                        tpl.maskGetAddrData(),
                        shape);
            break;
        case 3:
            normxcorrmw(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(), tpl.wghtTplGetAddrData(),
                        ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(), tpl.wghtRefGetAddrData(),
                        res.getAddrData(),
                        tpl.maskGetAddrData(),
                        shape);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cRegCorr<T>::computePCC(const cData3XCorr<T>& tpl,
                             const cData3XCorr<T>& ref,
                                   cData3XCorr<T>& res)
{
    const std::string   funcName("void cRegCorr<T>::computePCC("
                                    "const cData3XCorr<T>& tpl, "
                                    "const cData3XCorr<T>& ref, "
                                    "cData3XCorr<T>& res)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    switch (tpl.getDimension()) {
        case 1:
            pcorr(tpl.getAddrData(),
                  ref.getAddrData(),
                  res.getAddrData(),
                  tpl.getNrow());
            break;
        case 2:
            pcorr(tpl.getAddrData(),
                  ref.getAddrData(),
                  res.getAddrData(),
                  tpl.getNrow(), tpl.getNcol());
            break;
        case 3:
            pcorr(tpl.getAddrData(),
                  ref.getAddrData(),
                  res.getAddrData(),
                  tpl.getNrow(), tpl.getNcol(), tpl.getNsec());
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cRegCorr<T>::computeREG(const cData3XCorr<T>& tpl,
                             const cData3XCorr<T>& ref,
                                   cData3XCorr<T>& res,
                             T overlapRatio)
{
    const std::string   funcName("void cRegCorr<T>::computeREG("
                                    "const cData3XCorr<T>& tpl, "
                                    "const cData3XCorr<T>& ref, "
                                    "cData3XCorr<T>& res, "
                                    "T overlapRatio)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");
    require(tpl.getSize() == tpl.maskGetSize(),
            funcName + ": tpl and its mask have different sizes");
    require(ref.getSize() == ref.maskGetSize(),
            funcName + ": ref and its mask have different sizes");

    res.memReAlloc(ref.getSize() + tpl.getSize() - 1);
    res.maskGetRef()->memReAlloc(res.getSize());

    switch (tpl.getDimension()) {
        case 1:
            regnormxcorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                          ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                          tpl.maskGetAddrData(),
                          ref.maskGetAddrData(),
                          res.getAddrData(),
                          res.maskGetAddrData(),
                          overlapRatio);
            break;
        case 2:
            regnormxcorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                          ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                          tpl.maskGetAddrData(),
                          ref.maskGetAddrData(),
                          res.getAddrData(),
                          res.maskGetAddrData(),
                          overlapRatio);
            break;
        case 3:
            regnormxcorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                          ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                          tpl.maskGetAddrData(),
                          ref.maskGetAddrData(),
                          res.getAddrData(),
                          res.maskGetAddrData(),
                          overlapRatio);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

#ifdef __GEM_USE_CUDA__
template <typename T>
T cRegCorr<T>::computeXCC(const cuData3<T>& tpl, const cuData3<T>& ref)
{
    const std::string   funcName("T cRegCorr<T>::computeXCC("
                                    "const cuData3<T>& tpl, "
                                    "const cuData3<T>& ref)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    return cuda_xcorr_direct(tpl.getAddrData(),
                             ref.getAddrData(),
                             tpl.getNelement());
}
#endif

#ifdef __GEM_USE_CUDA__
template <typename T>
void cRegCorr<T>::computeXCC(const cuData3<T>& tpl,
                             const cuData3<T>& ref,
                                   cuData3<T>& res,
                             eXCorrRes shape)
{
    const std::string   funcName("void cRegCorr<T>::computeXCC("
                                    "const cuData3<T>& tpl, "
                                    "const cuData3<T>& ref, "
                                    "cuData3<T>& res, "
                                    "eXCorrRes shape)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() <= ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    switch (shape) {
        case XCORR_RES_FULL:
            res.memReAlloc(ref.getSize() + tpl.getSize() - cVector3<size_t>(1,1,1));
            break;
        case XCORR_RES_VALID:
            res.memReAlloc(ref.getSize() - tpl.getSize() + cVector3<size_t>(1,1,1));
            break;
        default:
            ERROR(funcName, "unsupported xcorr mode");
    }

    switch (tpl.getDimension()) {
        case 1:
            cuda_xcorr(tpl.getAddrData(), tpl.getNrow(),
                       ref.getAddrData(), ref.getNrow(),
                       res.getAddrData(),
                       shape);
            break;
        case 2:
            cuda_xcorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(),
                       ref.getAddrData(), ref.getNrow(), ref.getNcol(),
                       res.getAddrData(),
                       shape);
            break;
        case 3:
            cuda_xcorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                       ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                       res.getAddrData(),
                       shape);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}
#endif

#ifdef __GEM_USE_CUDA__
template <typename T>
T cRegCorr<T>::computeECC(const cuData3<T>& tpl,
                          const cuData3<T>& ref)
{
    const std::string   funcName("T cRegCorr<T>::computeECC("
                                    "const cuData3<T>& tpl, "
                                    "const cuData3<T>& ref)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    return cuda_ecorr_direct(tpl.getAddrData(),
                             ref.getAddrData(),
                             tpl.getNelement());

}
#endif

#ifdef __GEM_USE_CUDA__
template <typename T>
T cRegCorr<T>::computeECC(const cuData3<T>& tpl,
                          const cuData3<T>& ref,
                          const cuData3<T>& msk)
{
    const std::string   funcName("T cRegCorr<T>::computeECC("
                                    "const cuData3<T>& tpl, "
                                    "const cuData3<T>& ref, "
                                    "const cuData3<T>& msk)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");
    require(tpl.getSize() == msk.getSize(),
            funcName + ": tpl and msk have different sizes");

    return cuda_ecorrm_direct(tpl.getAddrData(),
                              ref.getAddrData(),
                              ref.getAddrData(),
                              tpl.getNelement());
}
#endif

#ifdef __GEM_USE_CUDA__
template <typename T>
void cRegCorr<T>::computeECC(const cuData3<T>& tpl,
                             const cuData3<T>& ref,
                                   cuData3<T>& res,
                             eXCorrRes shape)
{
    const std::string   funcName("void cRegCorr<T>::computeECC("
                                    "const cuData3<T>& tpl, "
                                    "const cuData3<T>& ref, "
                                    "cuData3<T>& res, "
                                    "eXCorrRes shape)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() <= ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    switch (shape) {
        case XCORR_RES_FULL:
            res.memReAlloc(ref.getSize() + tpl.getSize() - cVector3<size_t>(1,1,1));
            break;
        case XCORR_RES_VALID:
            res.memReAlloc(ref.getSize() - tpl.getSize() + cVector3<size_t>(1,1,1));
            break;
        default:
            ERROR(funcName, "unsupported xcorr mode");
    }

    switch (tpl.getDimension()) {
        case 1:
            cuda_ecorr(tpl.getAddrData(), tpl.getNrow(),
                       ref.getAddrData(), ref.getNrow(),
                       res.getAddrData(),
                       shape);
            break;
        case 2:
            cuda_ecorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(),
                       ref.getAddrData(), ref.getNrow(), ref.getNcol(),
                       res.getAddrData(),
                       shape);
            break;
        case 3:
            cuda_ecorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                       ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                       res.getAddrData(),
                       shape);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}
#endif

#ifdef __GEM_USE_CUDA__

template <typename T>
void cRegCorr<T>::computeECC(const cuData3<T>& tpl,
                             const cuData3<T>& ref,
                             const cuData3<T>& msk,
                                   cuData3<T>& res,
                             eXCorrRes shape)
{
    const std::string   funcName("void cRegCorr<T>::computeECC("
                                    "const cuData3<T>& tpl, "
                                    "const cuData3<T>& ref, "
                                    "const cuData3<T>& msk, "
                                    "cuData3<T>& res, "
                                    "eXCorrRes shape)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() <= ref.getSize(),
            funcName + ": tpl and ref have different sizes");
    require(tpl.getSize() == msk.getSize(),
            funcName + ": tpl and msk have different sizes");

    switch (shape) {
        case XCORR_RES_FULL:
            res.memReAlloc(ref.getSize() + tpl.getSize() - cVector3<size_t>(1,1,1));
            break;
        case XCORR_RES_VALID:
            res.memReAlloc(ref.getSize() - tpl.getSize() + cVector3<size_t>(1,1,1));
            break;
        default:
            ERROR(funcName, "unsupported xcorr mode");
    }

    switch (tpl.getDimension()) {
        case 1:
            cuda_ecorrm(tpl.getAddrData(), tpl.getNrow(),
                        ref.getAddrData(), ref.getNrow(),
                        res.getAddrData(),
                        msk.getAddrData(),
                        shape);
            break;
        case 2:
            cuda_ecorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(),
                        ref.getAddrData(), ref.getNrow(), ref.getNcol(),
                        res.getAddrData(),
                        msk.getAddrData(),
                        shape);
            break;
        case 3:
            cuda_ecorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                        ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                        res.getAddrData(),
                        msk.getAddrData(),
                        shape);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}
#endif

#ifdef __GEM_USE_CUDA__
template <typename T>
T cRegCorr<T>::computeNCC(const cuData3<T>& tpl,
                          const cuData3<T>& ref)
{
    const std::string   funcName("T cRegCorr<T>::computeNCC("
                                    "const cuData3<T>& tpl, "
                                    "const cuData3<T>& ref)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    return cuda_normxcorr_direct(tpl.getAddrData(),
                                 ref.getAddrData(),
                                 tpl.getNelement());
}
#endif

#ifdef __GEM_USE_CUDA__
template <typename T>
T cRegCorr<T>::computeNCC(const cuData3<T>& tpl,
                          const cuData3<T>& ref,
                          const cuData3<T>& msk)
{
    const std::string   funcName("T cRegCorr<T>::computeNCC("
                                    "const cuData3<T>& tpl, "
                                    "const cuData3<T>& ref, "
                                    "const cuData3<T>& msk)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() == ref.getSize(),
            funcName + ": tpl and ref have different sizes");
    require(tpl.getSize() == msk.getSize(),
            funcName + ": tpl and msk have different sizes");

    return cuda_normxcorrm_direct(tpl.getAddrData(),
                                  ref.getAddrData(),
                                  msk.getAddrData(),
                                  tpl.getNelement());
}
#endif

#ifdef __GEM_USE_CUDA__
template <typename T>
void cRegCorr<T>::computeNCC(const cuData3<T>& tpl,
                             const cuData3<T>& ref,
                                   cuData3<T>& res,
                             eXCorrRes shape)
{
    const std::string   funcName("void cRegCorr<T>::computeNCC("
                                    "const cuData3<T>& tpl, "
                                    "const cuData3<T>& ref, "
                                    "cuData3<T>& res, "
                                    "eXCorrRes shape)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() <= ref.getSize(),
            funcName + ": tpl and ref have different sizes");

    switch (shape) {
        case XCORR_RES_FULL:
            res.memReAlloc(ref.getSize() + tpl.getSize() - cVector3<size_t>(1,1,1));
            break;
        case XCORR_RES_VALID:
            res.memReAlloc(ref.getSize() - tpl.getSize() + cVector3<size_t>(1,1,1));
            break;
        default:
            ERROR(funcName, "unsupported xcorr mode");
    }

    switch (tpl.getDimension()) {
        case 1:
            cuda_normxcorr(tpl.getAddrData(), tpl.getNrow(),
                           ref.getAddrData(), ref.getNrow(),
                           res.getAddrData(),
                           shape);
            break;
        case 2:
            cuda_normxcorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(),
                           ref.getAddrData(), ref.getNrow(), ref.getNcol(),
                           res.getAddrData(),
                           shape);
            break;
        case 3:
            cuda_normxcorr(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                           ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                           res.getAddrData(),
                           shape);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}
#endif

#ifdef __GEM_USE_CUDA__
template <typename T>
void cRegCorr<T>::computeNCC(const cuData3<T>& tpl,
                             const cuData3<T>& ref,
                             const cuData3<T>& msk,
                                   cuData3<T>& res,
                             eXCorrRes shape)
{
    const std::string   funcName("void cRegCorr<T>::computeNCC("
                                    "const cuData3<T>& tpl, "
                                    "const cuData3<T>& ref, "
                                    "const cuData3<T>& msk, "
                                    "cuData3<T>& res, "
                                    "eXCorrRes shape)");

    tpl.requireNonEmpty(funcName);
    require(tpl.getSize() <= ref.getSize(),
            funcName + ": tpl and ref have different sizes");
    require(tpl.getSize() == msk.getSize(),
            funcName + ": tpl and msk have different sizes");

    switch (shape) {
        case XCORR_RES_FULL:
            res.memReAlloc(ref.getSize() + tpl.getSize() - cVector3<size_t>(1,1,1));
            break;
        case XCORR_RES_VALID:
            res.memReAlloc(ref.getSize() - tpl.getSize() + cVector3<size_t>(1,1,1));
            break;
        default:
            ERROR(funcName, "unsupported xcorr mode");
    }

    switch (tpl.getDimension()) {
        case 1:
            cuda_normxcorrm(tpl.getAddrData(), tpl.getNrow(),
                            ref.getAddrData(), ref.getNrow(),
                            res.getAddrData(),
                            msk.getAddrData(),
                            shape);
            break;
        case 2:
            cuda_normxcorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(),
                            ref.getAddrData(), ref.getNrow(), ref.getNcol(),
                            res.getAddrData(),
                            msk.getAddrData(),
                            shape);
            break;
        case 3:
            cuda_normxcorrm(tpl.getAddrData(), tpl.getNrow(), tpl.getNcol(), tpl.getNsec(),
                            ref.getAddrData(), ref.getNrow(), ref.getNcol(), ref.getNsec(),
                            res.getAddrData(),
                            msk.getAddrData(),
                            shape);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}
#endif

// instantiation
template class cRegCorr<float >;
template class cRegCorr<double>;

} // namespace gem

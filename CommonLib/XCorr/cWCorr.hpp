/***********************************************************************
 *  File:       cWCorr.hpp
 *
 *  Purpose:    Header file for WCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CWCORR_HPP__
#define __GEM_CWCORR_HPP__

#include "array.hpp"
#include "transform.hpp"
#include "fftw_lock.hpp"
#include "xcorr.hpp"

namespace gem {

class cWCorrSingle
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    float           *tplData,
                    *mskData,
                    *tplWght,
                    *refWght;

    float           *tplDataPadded,
                    *refDataPadded,
                    *refWghtPadded,
                    *resDataPadded;

    fftwf_complex   *tplDataPaddedFFT,
                    *refDataPaddedFFT,
                    *refDataPaddedFFT2,
                    *refWghtPaddedFFT,
                    *resDataPaddedFFT;

    float           *fftTmp1,
                    *fftTmp2;

    fftwf_plan      fftPlanTpl,
                    fftPlanRef,
                    fftPlanRef2,
                    fftPlanMsk,
                    fftPlanInv1,
                    fftPlanInv2,
                    fftPlanInv3;

    float           *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    bool            bMskSingle, bMskComputed;

    // 3D
    float           *tplDataRotInit,
                    *mskDataRotInit,
                    *tplWghtRotInit,
                    *refWghtRotInit;

public:
     cWCorrSingle(size_t refNrowInput,
                  size_t refNcolInput,
                  size_t refNsecInput = 1);
    ~cWCorrSingle();

    void setMskSingle(void) {bMskSingle = true;}

    void memAlloc(void);
    void memFree(void);

    void prepareRef(float *refDataInput);
    void setSizeTpl(size_t tplNrowInput,
                    size_t tplNcolInput,
                    size_t tplNsecInput = 1);
    void copyTplMskWght(float *tplDataInput, float *mskDataInput,
                        float *tplWghtInput, float *refWghtInput);
    void normalizeTpl();
    void prepareTplRefWght(void);

    void computeCorr(void);

    float  getResultMax(void);
    float* getResultRef(void);
    void   getResult(float *resDataPaddedOut);

    void mergeResult(size_t indx, eXCorrMerge bAbs);
    void mergeResultGlobal(float  *resDataAbsMaxPaddedGlobal,
                           size_t *resDataMaxIndPaddedGlobal,
                           eXCorrMerge bAbs);

    void printSize(void) {
        PRINT(tplNrow);    PRINT(tplNcol);    PRINT(tplNsec);
        PRINT(refNrow);    PRINT(refNcol);    PRINT(refNsec);
        PRINT(tplSize);    PRINT(refSize);    PRINT(fftSize);
    }

    float* getAddrTpl(void)     { return tplData; }
    float* getAddrMsk(void)     { return mskData; }
    float* getAddrTplWght(void) { return tplWght; }
    float* getAddrRefWght(void) { return refWght; }

    // 3D
    void copyTplMskWghtRot3D(float *tplDataInput, float *mskDataInput,
                             float *tplWghtInput, float *refWghtInput);

    void tplCreateByRot3D    (float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
    void mskCreateByRot3D    (float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
    void wghtTplCreateByRot3D(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
    void wghtRefCreateByRot3D(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
};

class cWCorrDouble
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    double          *tplData,
                    *mskData,
                    *tplWght,
                    *refWght;

    double          *tplDataPadded,
                    *refDataPadded,
                    *refWghtPadded,
                    *resDataPadded;

    fftw_complex    *tplDataPaddedFFT,
                    *refDataPaddedFFT,
                    *refDataPaddedFFT2,
                    *refWghtPaddedFFT,
                    *resDataPaddedFFT;

    double          *fftTmp1,
                    *fftTmp2;

    fftw_plan       fftPlanTpl,
                    fftPlanRef,
                    fftPlanRef2,
                    fftPlanMsk,
                    fftPlanInv1,
                    fftPlanInv2,
                    fftPlanInv3;

    double          *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    bool            bMskSingle, bMskComputed;

    // 3D
    double          *tplDataRotInit,
                    *mskDataRotInit,
                    *tplWghtRotInit,
                    *refWghtRotInit;

public:
     cWCorrDouble(size_t refNrowInput,
                  size_t refNcolInput,
                  size_t refNsecInput = 1);
    ~cWCorrDouble();

    void setMskSingle(void) {bMskSingle = true;}

    void memAlloc(void);
    void memFree(void);

    void prepareRef(double *refDataInput);
    void setSizeTpl(size_t tplNrowInput,
                    size_t tplNcolInput,
                    size_t tplNsecInput = 1);
    void copyTplMskWght(double *tplDataInput, double *mskDataInput,
                        double *tplWghtInput, double *refWghtInput);
    void normalizeTpl();
    void prepareTplRefWght(void);

    void computeCorr(void);

    double  getResultMax(void);
    double* getResultRef(void);
    void    getResult(double *resDataPaddedOut);

    void mergeResult(size_t indx, eXCorrMerge bAbs);
    void mergeResultGlobal(double *resDataAbsMaxPaddedGlobal,
                           size_t *resDataMaxIndPaddedGlobal,
                           eXCorrMerge bAbs);

    void printSize(void) {
        PRINT(tplNrow);    PRINT(tplNcol);    PRINT(tplNsec);
        PRINT(refNrow);    PRINT(refNcol);    PRINT(refNsec);
        PRINT(tplSize);    PRINT(refSize);    PRINT(fftSize);
    }

    double* getAddrTpl(void)     { return tplData; }
    double* getAddrMsk(void)     { return mskData; }
    double* getAddrTplWght(void) { return tplWght; }
    double* getAddrRefWght(void) { return refWght; }

    // 3D
    void copyTplMskWghtRot3D(double *tplDataInput, double *mskDataInput,
                             double *tplWghtInput, double *refWghtInput);

    void tplCreateByRot3D    (double alpha, double beta, double gamma, eInter inter = INTER_LINEAR);
    void mskCreateByRot3D    (double alpha, double beta, double gamma, eInter inter = INTER_LINEAR);
    void wghtTplCreateByRot3D(double alpha, double beta, double gamma, eInter inter = INTER_LINEAR);
    void wghtRefCreateByRot3D(double alpha, double beta, double gamma, eInter inter = INTER_LINEAR);
};

} // namespace gem

#endif

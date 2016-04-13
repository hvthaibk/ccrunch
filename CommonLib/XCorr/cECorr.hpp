/***********************************************************************
 *  File:       cECorr.hpp
 *
 *  Purpose:    Header file for ECC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CECORR_HPP__
#define __GEM_CECORR_HPP__

#include "array.hpp"
#include "transform.hpp"
#include "fftw_lock.hpp"
#include "xcorr.hpp"

namespace gem {

class cECorrSingle
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    float           *tplData,
                    *mskData;

    float           *tplDataPadded,
                    *refDataPadded,
                    *mskDataPadded,
                    *resDataPadded;

    fftwf_complex   *tplDataPaddedFFT,
                    *refDataPaddedFFT,
                    *refDataPaddedFFT2,
                    *mskDataPaddedFFT,
                    *resDataPaddedFFT;

    float           *fftTmp1;

    fftwf_plan      fftPlanTpl,
                    fftPlanRef,
                    fftPlanRef2,
                    fftPlanMsk,
                    fftPlanInv1,
                    fftPlanInv2;

    float           *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    bool            bMskSingle, bMskComputed;

    // 3D
    float           *tplDataRotInit,
                    *mskDataRotInit;

public:
     cECorrSingle(size_t refNrowInput,
                  size_t refNcolInput,
                  size_t refNsecInput = 1);
    ~cECorrSingle();

    void setMskSingle(void) {bMskSingle = true;}

    void memAlloc(void);
    void memFree(void);

    void prepareRef(float *refDataInput);
    void setSizeTpl(size_t tplNrowInput,
                    size_t tplNcolInput,
                    size_t tplNsecInput = 1);
    void copyTplMsk(float *tplDataInput, float *mskDataInput);
    void normalizeTpl(void);
    void prepareTplMsk(void);

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

    float* getAddrTpl(void) { return tplData; }
    float* getAddrMsk(void) { return mskData; }

    // 3D
    void copyTplMskRot3D(float *tplDataInput, float *mskDataInput);

    void tplCreateByRot3D(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
    void mskCreateByRot3D(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
};

class cECorrDouble
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    double          *tplData,
                    *mskData;

    double          *tplDataPadded,
                    *refDataPadded,
                    *mskDataPadded,
                    *resDataPadded;

    fftw_complex    *tplDataPaddedFFT,
                    *refDataPaddedFFT,
                    *refDataPaddedFFT2,
                    *mskDataPaddedFFT,
                    *resDataPaddedFFT;

    double          *fftTmp1;

    fftw_plan       fftPlanTpl,
                    fftPlanRef,
                    fftPlanRef2,
                    fftPlanMsk,
                    fftPlanInv1,
                    fftPlanInv2;

    double          *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    bool            bMskSingle, bMskComputed;

    // 3D
    double          *tplDataRotInit,
                    *mskDataRotInit;

public:
     cECorrDouble(size_t refNrowInput,
                  size_t refNcolInput,
                  size_t refNsecInput = 1);
    ~cECorrDouble();

    void setMskSingle(void) {bMskSingle = true;}

    void memAlloc(void);
    void memFree(void);

    void prepareRef(double *refDataInput);
    void setSizeTpl(size_t tplNrowInput,
                    size_t tplNcolInput,
                    size_t tplNsecInput = 1);
    void copyTplMsk(double *tplDataInput, double *mskDataInput);
    void normalizeTpl(void);
    void prepareTplMsk(void);

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

    double* getAddrTpl(void) { return tplData; }
    double* getAddrMsk(void) { return mskData; }

    // 3D
    void copyTplMskRot3D(double *tplDataInput, double *mskDataInput);

    void tplCreateByRot3D(double alpha, double beta, double gamma, eInter inter = INTER_LINEAR);
    void mskCreateByRot3D(double alpha, double beta, double gamma, eInter inter = INTER_LINEAR);
};

} // namespace gem

#endif

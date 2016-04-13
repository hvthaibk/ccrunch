/***********************************************************************
 *  File:       cXCorr.hpp
 *
 *  Purpose:    Header file for XCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CXCORR_HPP__
#define __GEM_CXCORR_HPP__

#include "array.hpp"
#include "transform.hpp"
#include "fftw_lock.hpp"
#include "xcorr.hpp"

namespace gem {

class cXCorrSingle
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    float           *tplData;

    float           *tplDataPadded,
                    *refDataPadded,
                    *resDataPadded;

    fftwf_complex   *tplDataPaddedFFT,
                    *refDataPaddedFFT,
                    *resDataPaddedFFT;

    fftwf_plan      fftPlanTpl,
                    fftPlanRef,
                    fftPlanInv;

    float           *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    // 3D
    float           *tplDataRotInit;

public:
     cXCorrSingle(size_t refNrowInput,
                  size_t refNcolInput,
                  size_t refNsecInput = 1);
    ~cXCorrSingle();

    void memAlloc(void);
    void memFree(void);

    void prepareRef(float *refDataInput);
    void setSizeTpl(size_t tplNrowInput,
                    size_t tplNcolInput,
                    size_t tplNsecInput = 1);
    void copyTpl(float *tplDataInput);
    void prepareTpl(void);

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

    // 3D
    void copyTplRot3D(float *tplDataInput);

    void tplCreateByRot3D(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
};

class cXCorrDouble
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    double          *tplData;

    double          *tplDataPadded,
                    *refDataPadded,
                    *resDataPadded;

    fftw_complex    *tplDataPaddedFFT,
                    *refDataPaddedFFT,
                    *resDataPaddedFFT;

    fftw_plan       fftPlanTpl,
                    fftPlanRef,
                    fftPlanInv;

    double          *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    // 3D
    double          *tplDataRotInit;

public:
     cXCorrDouble(size_t refNrowInput,
                  size_t refNcolInput,
                  size_t refNsecInput = 1);
    ~cXCorrDouble();

    void memAlloc(void);
    void memFree(void);

    void prepareRef(double *refDataInput);
    void setSizeTpl(size_t tplNrowInput,
                    size_t tplNcolInput,
                    size_t tplNsecInput = 1);
    void copyTpl(double *tplDataInput);
    void prepareTpl(void);

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

    // 3D
    void copyTplRot3D(double *tplDataInput);

    void tplCreateByRot3D(double alpha, double beta, double gamma, eInter inter = INTER_LINEAR);
};

} // namespace gem

#endif

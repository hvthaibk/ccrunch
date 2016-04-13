/***********************************************************************
 *  File:       cuECorr.cuh
 *
 *  Purpose:    Header file for ECC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUECORR_CUH__
#define __GEM_CUECORR_CUH__

#include "array.cuh"
#include "transform.cuh"
#include "xcorr.cuh"
#include <cufft.h>

namespace gem {

class cuECorrSingle
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    float           *dev_tplData,
                    *dev_mskData;

    float           *dev_tplDataPadded,
                    *dev_refDataPadded,
                    *dev_mskDataPadded,
                    *dev_resDataPadded;

    cuFloatComplex  *dev_tplDataPaddedFFT,
                    *dev_refDataPaddedFFT,
                    *dev_refDataPaddedFFT2,
                    *dev_mskDataPaddedFFT,
                    *dev_resDataPaddedFFT;

    float           *dev_fftTmp1;

    cufftHandle     fftPlanFwd,
                    fftPlanInv;

    float           *dev_resDataAbsMaxPadded;
    size_t          *dev_resDataMaxIndPadded;

    float           *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    bool            bMskSingle, bMskComputed;

    // 3D
    bool            bRotTexture;
    float           *dev_tplDataRotInit,
                    *dev_mskDataRotInit;
    cudaArray       *dev_tplDataRotInit_tex3D,
                    *dev_mskDataRotInit_tex3D;

public:
     cuECorrSingle(size_t refNrowInput,
                   size_t refNcolInput,
                   size_t refNsecInput = 1);
    ~cuECorrSingle();

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

    float* getAddrTpl(void) { return dev_tplData; }
    float* getAddrMsk(void) { return dev_mskData; }

    // 3D
    void setRotTexture(bool bRotTextureInput) { bRotTexture = bRotTextureInput; }
    bool getRotTexture(void                 ) { return bRotTexture;             }

    void copyTplMskRot3D(float *tplDataInput, float *mskDataInput);

    void tplCreateByRot3D(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
    void mskCreateByRot3D(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
};

class cuECorrDouble
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    double          *dev_tplData,
                    *dev_mskData;

    double          *dev_tplDataPadded,
                    *dev_refDataPadded,
                    *dev_mskDataPadded,
                    *dev_resDataPadded;

    cuDoubleComplex *dev_tplDataPaddedFFT,
                    *dev_refDataPaddedFFT,
                    *dev_refDataPaddedFFT2,
                    *dev_mskDataPaddedFFT,
                    *dev_resDataPaddedFFT;

    double          *dev_fftTmp1;

    cufftHandle     fftPlanFwd,
                    fftPlanInv;

    double          *dev_resDataAbsMaxPadded;
    size_t          *dev_resDataMaxIndPadded;

    double          *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    bool            bMskSingle, bMskComputed;

    // 3D
    double          *dev_tplDataRotInit,
                    *dev_mskDataRotInit;

public:
     cuECorrDouble(size_t refNrowInput,
                   size_t refNcolInput,
                   size_t refNsecInput = 1);
    ~cuECorrDouble();

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

    double* getAddrTpl(void) { return dev_tplData; }
    double* getAddrMsk(void) { return dev_mskData; }

    // 3D
    void copyTplMskRot3D(double *tplDataInput, double *mskDataInput);

    void tplCreateByRot3D(double alpha, double beta, double gamma, eInter inter = INTER_LINEAR);
    void mskCreateByRot3D(double alpha, double beta, double gamma, eInter inter = INTER_LINEAR);
};

} // namespace gem

#endif

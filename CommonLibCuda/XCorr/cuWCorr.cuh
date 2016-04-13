/***********************************************************************
 *  File:       cuWCorr.cuh
 *
 *  Purpose:    Header file for WCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUWCORR_CUH__
#define __GEM_CUWCORR_CUH__

#include "array.cuh"
#include "transform.cuh"
#include "xcorr.cuh"
#include <cufft.h>

namespace gem {

class cuWCorrSingle
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    float           *dev_tplData,
                    *dev_mskData,
                    *dev_tplWght,
                    *dev_refWght;

    float           *dev_tplDataPadded,
                    *dev_refDataPadded,
                    *dev_refWghtPadded,
                    *dev_resDataPadded;

    cuFloatComplex  *dev_tplDataPaddedFFT,
                    *dev_refDataPaddedFFT,
                    *dev_refDataPaddedFFT2,
                    *dev_refWghtPaddedFFT,
                    *dev_resDataPaddedFFT;

    float           *dev_fftTmp1,
                    *dev_fftTmp2;

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
                    *dev_mskDataRotInit,
                    *dev_tplWghtRotInit,
                    *dev_refWghtRotInit;
    cudaArray       *dev_tplDataRotInit_tex3D,
                    *dev_mskDataRotInit_tex3D,
                    *dev_tplWghtRotInit_tex3D,
                    *dev_refWghtRotInit_tex3D;

public:
     cuWCorrSingle(size_t refNrowInput,
                   size_t refNcolInput,
                   size_t refNsecInput = 1);
    ~cuWCorrSingle();

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
    void prepareTplRefWght();

    void computeCorr();

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

    float* getAddrTpl()     { return dev_tplData; }
    float* getAddrMsk()     { return dev_mskData; }
    float* getAddrTplWght() { return dev_tplWght; }
    float* getAddrRefWght() { return dev_refWght; }

    // 3D
    void setRotTexture(bool bRotTextureInput) { bRotTexture = bRotTextureInput; }
    bool getRotTexture(void                 ) { return bRotTexture;             }

    void copyTplMskWghtRot3D(float *tplDataInput, float *mskDataInput,
                             float *tplWghtInput, float *refWghtInput);

    void tplCreateByRot3D    (float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
    void mskCreateByRot3D    (float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
    void wghtTplCreateByRot3D(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
    void wghtRefCreateByRot3D(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
};

class cuWCorrDouble
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    double          *dev_tplData,
                    *dev_mskData,
                    *dev_tplWght,
                    *dev_refWght;

    double          *dev_tplDataPadded,
                    *dev_refDataPadded,
                    *dev_refWghtPadded,
                    *dev_resDataPadded;

    cuDoubleComplex *dev_tplDataPaddedFFT,
                    *dev_refDataPaddedFFT,
                    *dev_refDataPaddedFFT2,
                    *dev_refWghtPaddedFFT,
                    *dev_resDataPaddedFFT;

    double          *dev_fftTmp1,
                    *dev_fftTmp2;

    cufftHandle     fftPlanFwd,
                    fftPlanInv;

    double          *dev_resDataAbsMaxPadded;
    size_t          *dev_resDataMaxIndPadded;

    double          *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    bool            bMskSingle, bMskComputed;

    // 3D
    double          *dev_tplDataRotInit,
                    *dev_mskDataRotInit,
                    *dev_tplWghtRotInit,
                    *dev_refWghtRotInit;

public:
     cuWCorrDouble(size_t refNrowInput,
                   size_t refNcolInput,
                   size_t refNsecInput = 1);
    ~cuWCorrDouble();

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
    void prepareTplRefWght();

    void computeCorr();

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

    double* getAddrTpl()     { return dev_tplData; }
    double* getAddrMsk()     { return dev_mskData; }
    double* getAddrTplWght() { return dev_tplWght; }
    double* getAddrRefWght() { return dev_refWght; }

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

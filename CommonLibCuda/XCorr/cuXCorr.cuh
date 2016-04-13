/***********************************************************************
 *  File:       cuXCorr.cuh
 *
 *  Purpose:    Header file for XCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUXCORR_CUH__
#define __GEM_CUXCORR_CUH__

#include "array.cuh"
#include "transform.cuh"
#include "xcorr.cuh"
#include <cufft.h>

namespace gem {

class cuXCorrSingle
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    float           *dev_tplData;

    float           *dev_tplDataPadded,
                    *dev_refDataPadded,
                    *dev_resDataPadded;

    cuFloatComplex  *dev_tplDataPaddedFFT,
                    *dev_refDataPaddedFFT,
                    *dev_resDataPaddedFFT;

    cufftHandle     fftPlanFwd,
                    fftPlanInv;

    float           *dev_resDataAbsMaxPadded;
    size_t          *dev_resDataMaxIndPadded;

    float           *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    // 3D
    bool            bRotTexture;
    float           *dev_tplDataRotInit;
    cudaArray       *dev_tplDataRotInit_tex3D;

public:
     cuXCorrSingle(size_t refNrowInput,
                   size_t refNcolInput,
                   size_t refNsecInput = 1);
    ~cuXCorrSingle();

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

    float* getAddrTpl(void) { return dev_tplData; }

    // 3D
    void setRotTexture(bool bRotTextureInput) { bRotTexture = bRotTextureInput; }
    bool getRotTexture(void                 ) { return bRotTexture;             }

    void copyTplRot3D(float *tplDataInput);

    void tplCreateByRot3D(float alpha, float beta, float gamma, eInter inter = INTER_LINEAR);
};

class cuXCorrDouble
{
private:
    size_t          tplNrow, tplNcol, tplNsec,
                    refNrow, refNcol, refNsec;
    size_t          tplSize, refSize, fftSize;

    double          *dev_tplData;

    double          *dev_tplDataPadded,
                    *dev_refDataPadded,
                    *dev_resDataPadded;

    cuDoubleComplex *dev_tplDataPaddedFFT,
                    *dev_refDataPaddedFFT,
                    *dev_resDataPaddedFFT;

    cufftHandle     fftPlanFwd,
                    fftPlanInv;

    double          *dev_resDataAbsMaxPadded;
    size_t          *dev_resDataMaxIndPadded;

    double          *resDataAbsMaxPadded;
    size_t          *resDataMaxIndPadded;

    // 3D
    double          *dev_tplDataRotInit;

public:
     cuXCorrDouble(size_t refNrowInput,
                   size_t refNcolInput,
                   size_t refNsecInput = 1);
    ~cuXCorrDouble();

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

    double* getAddrTpl(void) { return dev_tplData; }

    // 3D
    void copyTplRot3D(double *tplDataInput);

    void tplCreateByRot3D(double alpha, double beta, double gamma, eInter inter = INTER_LINEAR);
};

} // namespace gem

#endif

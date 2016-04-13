/***********************************************************************
 *  File:       cFitterRefineGSL.cpp
 *
 *  Purpose:    Implementation of cFitter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cFitter.hpp"

#include "gsl/gsl_multimin.h"

cData3<double>      powRef, powTpl,    powMsk,    powWghtRef,    powWghtTpl,
                            powTplTmp, powMskTmp, powWghtRefTmp, powWghtTplTmp;
#ifdef __GEM_USE_CUDA__
cuData3<double>     powRefCUDA, powTplCUDA,    powMskCUDA,    powWghtRefCUDA,    powWghtTplCUDA,
                                powTplCUDATmp, powMskCUDATmp, powWghtRefCUDATmp, powWghtTplCUDATmp;
#endif

double gslFuncXCorr(const gsl_vector *v, void *)
{
    double x[6];

    x[0] = gsl_vector_get (v, 0);
    x[1] = gsl_vector_get (v, 1);
    x[2] = gsl_vector_get (v, 2);
    x[3] = gsl_vector_get (v, 3);
    x[4] = gsl_vector_get (v, 4);
    x[5] = gsl_vector_get (v, 5);

    transform_rotate_translate(powTpl.getAddrData(), powTplTmp.getAddrData(),
                               powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                               x[0], x[1], x[2],
                               x[3], x[4], x[5]);

    return -xcorr_direct(powTplTmp.getAddrData(),
                         powRef   .getAddrData(),
                         powTplTmp.getNelement());
}

#ifdef __GEM_USE_CUDA__
double gslFuncXCorrCUDA(const gsl_vector *v, void *)
#else
double gslFuncXCorrCUDA(const gsl_vector * , void *)
#endif
{
#ifdef __GEM_USE_CUDA__
    double x[6];

    x[0] = gsl_vector_get (v, 0);
    x[1] = gsl_vector_get (v, 1);
    x[2] = gsl_vector_get (v, 2);
    x[3] = gsl_vector_get (v, 3);
    x[4] = gsl_vector_get (v, 4);
    x[5] = gsl_vector_get (v, 5);

    cuda_transform_rotate_translate(powTplCUDA.getAddrData(), powTplCUDATmp.getAddrData(),
                                    powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                                    x[0], x[1], x[2],
                                    x[3], x[4], x[5]);

    return -cuda_xcorr_direct(powTplCUDATmp.getAddrData(),
                              powRefCUDA   .getAddrData(),
                              powTplCUDATmp.getNelement());
#else
    std::cout << "gslFuncXCorrCUDA(): "
              << "this function is not supported without CUDA"
              << std::endl;
    return 0;
#endif
}

double gslFuncECorr(const gsl_vector *v, void *)
{
    double x[6];

    x[0] = gsl_vector_get (v, 0);
    x[1] = gsl_vector_get (v, 1);
    x[2] = gsl_vector_get (v, 2);
    x[3] = gsl_vector_get (v, 3);
    x[4] = gsl_vector_get (v, 4);
    x[5] = gsl_vector_get (v, 5);

    transform_rotate_translate(powTpl.getAddrData(), powTplTmp.getAddrData(),
                               powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                               x[0], x[1], x[2],
                               x[3], x[4], x[5]);

    transform_rotate_translate(powMsk.getAddrData(), powMskTmp.getAddrData(),
                               powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                               x[0], x[1], x[2],
                               x[3], x[4], x[5]);

    cTransThreshold<double>().threshold(powMskTmp, 0.5f);

    return -ecorrm_direct(powTplTmp.getAddrData(),
                          powRef   .getAddrData(),
                          powMskTmp.getAddrData(),
                          powTplTmp.getNelement());
}

#ifdef __GEM_USE_CUDA__
double gslFuncECorrCUDA(const gsl_vector *v, void *)
#else
double gslFuncECorrCUDA(const gsl_vector * , void *)
#endif
{
#ifdef __GEM_USE_CUDA__
    double x[6];

    x[0] = gsl_vector_get (v, 0);
    x[1] = gsl_vector_get (v, 1);
    x[2] = gsl_vector_get (v, 2);
    x[3] = gsl_vector_get (v, 3);
    x[4] = gsl_vector_get (v, 4);
    x[5] = gsl_vector_get (v, 5);

    cuda_transform_rotate_translate(powTplCUDA.getAddrData(), powTplCUDATmp.getAddrData(),
                                    powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                                    x[0], x[1], x[2],
                                    x[3], x[4], x[5]);

    cuda_transform_rotate_translate(powMskCUDA.getAddrData(), powMskCUDATmp.getAddrData(),
                                    powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                                    x[0], x[1], x[2],
                                    x[3], x[4], x[5]);

    cuTransThreshold<double>().threshold(powMskCUDATmp, 0.5f);

    return -cuda_ecorrm_direct(powTplCUDATmp.getAddrData(),
                               powRefCUDA   .getAddrData(),
                               powMskCUDATmp.getAddrData(),
                               powTplCUDATmp.getNelement());
#else
    std::cout << "gslFuncECorrCUDA(): "
              << "this function is not supported without CUDA"
              << std::endl;
    return 0;
#endif
}

double gslFuncNCorr(const gsl_vector *v, void *)
{
    double x[6];

    x[0] = gsl_vector_get (v, 0);
    x[1] = gsl_vector_get (v, 1);
    x[2] = gsl_vector_get (v, 2);
    x[3] = gsl_vector_get (v, 3);
    x[4] = gsl_vector_get (v, 4);
    x[5] = gsl_vector_get (v, 5);

    transform_rotate_translate(powTpl.getAddrData(), powTplTmp.getAddrData(),
                               powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                               x[0], x[1], x[2],
                               x[3], x[4], x[5]);

    transform_rotate_translate(powMsk.getAddrData(), powMskTmp.getAddrData(),
                               powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                               x[0], x[1], x[2],
                               x[3], x[4], x[5]);

    cTransThreshold<double>().threshold(powMskTmp, 0.5f);

    return -normxcorrm_direct(powTplTmp.getAddrData(),
                              powRef   .getAddrData(),
                              powMskTmp.getAddrData(),
                              powTplTmp.getNelement());
}

#ifdef __GEM_USE_CUDA__
double gslFuncNCorrCUDA(const gsl_vector *v, void *)
#else
double gslFuncNCorrCUDA(const gsl_vector * , void *)
#endif
{
#ifdef __GEM_USE_CUDA__
    double x[6];

    x[0] = gsl_vector_get (v, 0);
    x[1] = gsl_vector_get (v, 1);
    x[2] = gsl_vector_get (v, 2);
    x[3] = gsl_vector_get (v, 3);
    x[4] = gsl_vector_get (v, 4);
    x[5] = gsl_vector_get (v, 5);

    cuda_transform_rotate_translate(powTplCUDA.getAddrData(), powTplCUDATmp.getAddrData(),
                                    powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                                    x[0], x[1], x[2],
                                    x[3], x[4], x[5]);

    cuda_transform_rotate_translate(powMskCUDA.getAddrData(), powMskCUDATmp.getAddrData(),
                                    powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                                    x[0], x[1], x[2],
                                    x[3], x[4], x[5]);

    cuTransThreshold<double>().threshold(powMskCUDATmp, 0.5f);

    return -cuda_normxcorrm_direct(powTplCUDATmp.getAddrData(),
                                   powRefCUDA   .getAddrData(),
                                   powMskCUDATmp.getAddrData(),
                                   powTplCUDATmp.getNelement());
#else
    std::cout << "gslFuncNCorrCUDA(): "
              << "this function is not supported without CUDA"
              << std::endl;
    return 0;
#endif
}

double gslFuncWCorr(const gsl_vector *v, void *)
{
    double x[6];

    x[0] = gsl_vector_get (v, 0);
    x[1] = gsl_vector_get (v, 1);
    x[2] = gsl_vector_get (v, 2);
    x[3] = gsl_vector_get (v, 3);
    x[4] = gsl_vector_get (v, 4);
    x[5] = gsl_vector_get (v, 5);

    transform_rotate_translate(powTpl.getAddrData(), powTplTmp.getAddrData(),
                               powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                               x[0], x[1], x[2],
                               x[3], x[4], x[5]);

    transform_rotate_translate(powMsk.getAddrData(), powMskTmp.getAddrData(),
                               powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                               x[0], x[1], x[2],
                               x[3], x[4], x[5]);

    cTransThreshold<double>().threshold(powMskTmp, 0.5f);

    transform_rotate_translate(powWghtTpl.getAddrData(), powWghtTplTmp.getAddrData(),
                               powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                               x[0], x[1], x[2],
                               x[3], x[4], x[5]);

    transform_rotate_translate(powWghtRef.getAddrData(), powWghtRefTmp.getAddrData(),
                               powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                               x[0], x[1], x[2],
                               x[3], x[4], x[5]);

    return -normxcorrmw_direct(powTplTmp.getAddrData(), powWghtTplTmp.getAddrData(),
                               powRef   .getAddrData(), powWghtRefTmp.getAddrData(),
                               powMskTmp.getAddrData(),
                               powTplTmp.getNelement());
}

#ifdef __GEM_USE_CUDA__
double gslFuncWCorrCUDA(const gsl_vector *v, void *)
#else
double gslFuncWCorrCUDA(const gsl_vector * , void *)
#endif
{
#ifdef __GEM_USE_CUDA__
    double x[6];

    x[0] = gsl_vector_get (v, 0);
    x[1] = gsl_vector_get (v, 1);
    x[2] = gsl_vector_get (v, 2);
    x[3] = gsl_vector_get (v, 3);
    x[4] = gsl_vector_get (v, 4);
    x[5] = gsl_vector_get (v, 5);

    cuda_transform_rotate_translate(powTplCUDA.getAddrData(), powTplCUDATmp.getAddrData(),
                                    powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                                    x[0], x[1], x[2],
                                    x[3], x[4], x[5]);

    cuda_transform_rotate_translate(powMskCUDA.getAddrData(), powMskCUDATmp.getAddrData(),
                                    powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                                    x[0], x[1], x[2],
                                    x[3], x[4], x[5]);

    cuTransThreshold<double>().threshold(powMskCUDATmp, 0.5f);

    cuda_transform_rotate_translate(powWghtTplCUDA.getAddrData(), powWghtTplCUDATmp.getAddrData(),
                                    powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                                    x[0], x[1], x[2],
                                    x[3], x[4], x[5]);

    cuda_transform_rotate_translate(powWghtRefCUDA.getAddrData(), powWghtRefCUDATmp.getAddrData(),
                                    powTpl.getNrow(), powTpl.getNcol(), powTpl.getNsec(),
                                    x[0], x[1], x[2],
                                    x[3], x[4], x[5]);

    return -cuda_normxcorrmw_direct(powTplCUDATmp.getAddrData(), powWghtTplCUDATmp.getAddrData(),
                                    powRefCUDA   .getAddrData(), powWghtRefCUDATmp.getAddrData(),
                                    powMskCUDATmp.getAddrData(),
                                    powTplCUDATmp.getNelement());
#else
    std::cout << "gslFuncWCorrCUDA(): "
              << "this function is not supported without CUDA"
              << std::endl;
    return 0;
#endif
}

double cFitter::refineGSL(cVector3<double>& transDouble, cVector3<double>& angleDouble)
{
    const std::string   funcName("double cFitter::refineGSL("
                                "cVector3<double>& transDouble, "
                                "cVector3<double>& angleDouble)");

    switch (mCorr) {
        case CC_MODE_WCORR:
        case CC_MODE_CCORR:
            powWghtRefTmp.memReAlloc(powMsk.getSize());
            powWghtTplTmp.memReAlloc(powMsk.getSize());
        case CC_MODE_XCORR:
        case CC_MODE_ECORR:
        case CC_MODE_NCORR:
        case CC_MODE_SCORR:
            powTplTmp.memReAlloc(powTpl.getSize());
            powMskTmp.memReAlloc(powMsk.getSize());
            break;
        default:
            ERROR(funcName, "unsupported correlation mode");
    }

#ifdef __GEM_USE_CUDA__
    if (nGPU > 0) {
        switch (mCorr) {
            case CC_MODE_WCORR:
            case CC_MODE_CCORR:
                powWghtRefCUDATmp.memReAlloc(powMsk.getSize());
                powWghtTplCUDATmp.memReAlloc(powMsk.getSize());
            case CC_MODE_XCORR:
            case CC_MODE_ECORR:
            case CC_MODE_NCORR:
            case CC_MODE_SCORR:
                powTplCUDATmp.memReAlloc(powTpl.getSize());
                powMskCUDATmp.memReAlloc(powMsk.getSize());
                break;
            default:
                ERROR(funcName, "unsupported correlation mode");
        }
    }
#endif

    double                              funcMin = 0;
    const gsl_multimin_fminimizer_type  *T = gsl_multimin_fminimizer_nmsimplex2rand;
          gsl_multimin_fminimizer       *s = NULL;
          gsl_multimin_function         minex_func;
          gsl_vector                    *ss, *x;

    // Starting point
    x = gsl_vector_alloc (6);
    gsl_vector_set (x, 0, angleDouble[0]);
    gsl_vector_set (x, 1, angleDouble[1]);
    gsl_vector_set (x, 2, angleDouble[2]);
    gsl_vector_set (x, 3, 0);
    gsl_vector_set (x, 4, 0);
    gsl_vector_set (x, 5, 0);

    // Set initial step sizes to 1
    ss = gsl_vector_alloc (6);
    gsl_vector_set_all (ss, 1.0);

    // Initialize method and iterate
    minex_func.n = 6;
    minex_func.params = NULL;

    switch (mCorr) {
        case CC_MODE_XCORR:
            if (nGPU > 0) minex_func.f = gslFuncXCorrCUDA;
            else          minex_func.f = gslFuncXCorr;

            tolRefine = tolRefine*powRef.getSum2();
            break;
        case CC_MODE_ECORR:
            if (nGPU > 0) minex_func.f = gslFuncECorrCUDA;
            else          minex_func.f = gslFuncECorr;
            break;
        case CC_MODE_NCORR:
            if (nGPU > 0) minex_func.f = gslFuncNCorrCUDA;
            else          minex_func.f = gslFuncNCorr;
            break;
        case CC_MODE_WCORR:
            if (nGPU > 0) minex_func.f = gslFuncWCorrCUDA;
            else          minex_func.f = gslFuncWCorr;
            break;
        case CC_MODE_SCORR:
            ERROR(funcName, "unsupported correlation mode");
            break;
        case CC_MODE_CCORR:
            ERROR(funcName, "unsupported correlation mode");
            break;
        default:
            ERROR(funcName, "unsupported correlation mode");
    }

    s = gsl_multimin_fminimizer_alloc (T, 6);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    size_t  iter = 0;
    int     status;
    double  size;

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);
        size   = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, tolRefine);
    }
    while (status == GSL_CONTINUE && iter < iterRefine);

    angleDouble[0] = gsl_vector_get (s->x, 0);
    angleDouble[1] = gsl_vector_get (s->x, 1);
    angleDouble[2] = gsl_vector_get (s->x, 2);
    transDouble[0] = gsl_vector_get (s->x, 3);
    transDouble[1] = gsl_vector_get (s->x, 4);
    transDouble[2] = gsl_vector_get (s->x, 5);
    funcMin = s->fval;

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    powTplTmp.memFree();
    powMskTmp.memFree();
    powWghtRefTmp.memFree();
    powWghtTplTmp.memFree();

#ifdef __GEM_USE_CUDA__
    if (nGPU > 0) {
        powTplCUDATmp.memFree();
        powMskCUDATmp.memFree();
        powWghtRefCUDATmp.memFree();
        powWghtTplCUDATmp.memFree();
    }
#endif

    return funcMin;
}

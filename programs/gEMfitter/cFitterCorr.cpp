/***********************************************************************
 *  File:       cFitterCorr.cpp
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

#include "pft_vertices.hpp"

static cQueueConcurrent<cData3XCorr<float>* >   queueAngle(30);
static cQueueConcurrent<std::string>            queueMsge;
static bool                                     streamData, streamMsge;
static cData3XCorr<float>                       *mapRefLocal = NULL, *mapTplLocal = NULL;
static cData3XCorr<float>                       *mapRefSurfaceLocal = NULL, *mapRefOuterLocal = NULL;
static cData3XCorr<float>                       *mapTplSurfaceLocal = NULL, *mapTplInnerLocal = NULL;


static bool                                     bRotTextureLocal = true;
static float                                    wDensityLocal = 0.5f, wSurfaceLocal = 0.5f, wPenaltyLocal = 0.5f;

static cData3XCorr<float>                       resDataAbsMaxPaddedGlobal;
static cArray3<size_t>                          resDataMaxIndPaddedGlobal;

#ifndef NDEBUG
void threadFuncXCorr(int threadID, int mpiID)
#else
void threadFuncXCorr(int , int )
#endif
{
    std::stringstream                   strStream;
    boost::posix_time::microseconds     waitingTime(1);
    cData3XCorr<float>                  *angleCorr = NULL;

    // cXCorr object & prepare ref
    cXCorrSingle    cXCorr(mapRefLocal->getNrow(),
                           mapRefLocal->getNcol(),
                           mapRefLocal->getNsec());

    cXCorr.memAlloc();
    cXCorr.prepareRef(mapRefLocal->getAddrData());
    cXCorr.setSizeTpl(mapTplLocal->getNrow(),
                      mapTplLocal->getNcol(),
                      mapTplLocal->getNsec());
    cXCorr.copyTplRot3D(mapTplLocal->getAddrData());

    while (streamData || !queueAngle.empty()) {
        if (queueAngle.wait_and_pop(angleCorr, waitingTime)) {
/*#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    rotation"   << std::setw(9) << angleCorr->getID()
                      << ",    threadCPU " << std::setw(6) << threadID
                      << ",    angle1 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(0)
                      << ",    angle2 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(1)
                      << ",    angle3 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(2);
            queueMsge.push(strStream.str());
#endif*/

            // prepare tpl and msk
            cXCorr.tplCreateByRot3D(angleCorr->getElement(0),
                                    angleCorr->getElement(1),
                                    angleCorr->getElement(2),
                                    INTER_LINEAR);
            cXCorr.prepareTpl();

            // compute SCC
            cXCorr.computeCorr();

            // merge result
            cXCorr.mergeResult(angleCorr->getID(), XCORR_MERGE_NEGATIVE);   // negative for surface fitting

            delete angleCorr;
        }
    }

    // merge result global
    cXCorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                             resDataMaxIndPaddedGlobal.getAddrData(),
                             XCORR_MERGE_NEGATIVE);

    cXCorr.memFree();
}

#ifndef NDEBUG
void threadFuncXCorrCUDA(int threadID, int mpiID)
#else
#ifdef __GEM_USE_CUDA__
void threadFuncXCorrCUDA(int threadID, int )
#else
void threadFuncXCorrCUDA(int , int )
#endif
#endif
{
#ifdef __GEM_USE_CUDA__
    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());

    std::stringstream                   strStream;
    boost::posix_time::microseconds     waitingTime(1);
    cData3XCorr<float>                  *angleCorr = NULL;

    // cXCorr object & prepare ref
    cuXCorrSingle     cuXCorr(mapRefLocal->getNrow(),
                              mapRefLocal->getNcol(),
                              mapRefLocal->getNsec());

    cuXCorr.memAlloc();
    cuXCorr.prepareRef(mapRefLocal->getAddrData());
    cuXCorr.setSizeTpl(mapTplLocal->getNrow(),
                       mapTplLocal->getNcol(),
                       mapTplLocal->getNsec());
    cuXCorr.setRotTexture(bRotTextureLocal);
    cuXCorr.copyTplRot3D(mapTplLocal->getAddrData());

    while (streamData || !queueAngle.empty()) {
        if (queueAngle.wait_and_pop(angleCorr, waitingTime)) {
/*#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    rotation"   << std::setw(9) << angleCorr->getID()
                      << ",    threadGPU " << std::setw(6) << threadID
                      << ",    angle1 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(0)
                      << ",    angle2 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(1)
                      << ",    angle3 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(2);
            queueMsge.push(strStream.str());
#endif*/

            // prepare tpl and msk
            cuXCorr.tplCreateByRot3D(angleCorr->getElement(0),
                                     angleCorr->getElement(1),
                                     angleCorr->getElement(2),
                                     INTER_LINEAR);
            cuXCorr.prepareTpl();

            // compute SCC
            cuXCorr.computeCorr();

            // merge result
            cuXCorr.mergeResult(angleCorr->getID(), XCORR_MERGE_NEGATIVE);  // negative for surface fitting

            delete angleCorr;
        }
    }

    // merge result global
    cuXCorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                              resDataMaxIndPaddedGlobal.getAddrData(),
                              XCORR_MERGE_NEGATIVE);

    cuXCorr.memFree();

    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());
#else
    std::cout << "threadFuncXCorrCUDA(): "
              << "this function is not supported without CUDA"
              << std::endl;
#endif
}

#ifndef NDEBUG
void threadFuncECorr(int threadID, int mpiID)
#else
void threadFuncECorr(int , int )
#endif
{
    std::stringstream                   strStream;
    boost::posix_time::microseconds     waitingTime(1);
    cData3XCorr<float>                  *angleCorr = NULL;

    // cECorr object & prepare ref
    cECorrSingle    cECorr(mapRefLocal->getNrow(),
                           mapRefLocal->getNcol(),
                           mapRefLocal->getNsec());

    cECorr.memAlloc();
    cECorr.prepareRef(mapRefLocal->getAddrData());
    cECorr.setSizeTpl(mapTplLocal->getNrow(),
                      mapTplLocal->getNcol(),
                      mapTplLocal->getNsec());
    cECorr.copyTplMskRot3D(mapTplLocal->getAddrData(),
                           mapTplLocal->maskGetAddrData());

    while (streamData || !queueAngle.empty()) {
        if (queueAngle.wait_and_pop(angleCorr, waitingTime)) {
/*#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    rotation"   << std::setw(9) << angleCorr->getID()
                      << ",    threadCPU " << std::setw(6) << threadID
                      << ",    angle1 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(0)
                      << ",    angle2 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(1)
                      << ",    angle3 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(2);
            queueMsge.push(strStream.str());
#endif*/

            // prepare tpl and msk
            cECorr.tplCreateByRot3D(angleCorr->getElement(0),
                                    angleCorr->getElement(1),
                                    angleCorr->getElement(2),
                                    INTER_LINEAR);
            cECorr.mskCreateByRot3D(angleCorr->getElement(0),
                                    angleCorr->getElement(1),
                                    angleCorr->getElement(2),
                                    INTER_LINEAR);
            cECorr.normalizeTpl();
            cECorr.prepareTplMsk();

            // compute XCC
            cECorr.computeCorr();

            // merge result
            cECorr.mergeResult(angleCorr->getID(), XCORR_MERGE_POSITIVE);

            delete angleCorr;
        }
    }

    // merge result global
    cECorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                             resDataMaxIndPaddedGlobal.getAddrData(),
                             XCORR_MERGE_POSITIVE);

    cECorr.memFree();
}

#ifndef NDEBUG
void threadFuncECorrCUDA(int threadID, int mpiID)
#else
#ifdef __GEM_USE_CUDA__
void threadFuncECorrCUDA(int threadID, int )
#else
void threadFuncECorrCUDA(int , int )
#endif
#endif
{
#ifdef __GEM_USE_CUDA__
    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());

    std::stringstream                   strStream;
    boost::posix_time::microseconds     waitingTime(1);
    cData3XCorr<float>                  *angleCorr = NULL;

    // cECorr object & prepare ref
    cuECorrSingle       cECorrCUDA(mapRefLocal->getNrow(),
                                   mapRefLocal->getNcol(),
                                   mapRefLocal->getNsec());

    cECorrCUDA.memAlloc();
    cECorrCUDA.prepareRef(mapRefLocal->getAddrData());
    cECorrCUDA.setSizeTpl(mapTplLocal->getNrow(),
                          mapTplLocal->getNcol(),
                          mapTplLocal->getNsec());
    cECorrCUDA.setRotTexture(bRotTextureLocal);
    cECorrCUDA.copyTplMskRot3D(mapTplLocal->getAddrData(),
                               mapTplLocal->maskGetAddrData());

    while (streamData || !queueAngle.empty()) {
        if (queueAngle.wait_and_pop(angleCorr, waitingTime)) {
/*#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    rotation"   << std::setw(9) << angleCorr->getID()
                      << ",    threadGPU " << std::setw(6) << threadID
                      << ",    angle1 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(0)
                      << ",    angle2 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(1)
                      << ",    angle3 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(2);
            queueMsge.push(strStream.str());
#endif*/

            // prepare tpl and msk
            cECorrCUDA.tplCreateByRot3D(angleCorr->getElement(0),
                                        angleCorr->getElement(1),
                                        angleCorr->getElement(2),
                                        INTER_LINEAR);
            cECorrCUDA.mskCreateByRot3D(angleCorr->getElement(0),
                                        angleCorr->getElement(1),
                                        angleCorr->getElement(2),
                                        INTER_LINEAR);
            cECorrCUDA.normalizeTpl();
            cECorrCUDA.prepareTplMsk();

            // compute XCC
            cECorrCUDA.computeCorr();

            // merge result
            cECorrCUDA.mergeResult(angleCorr->getID(), XCORR_MERGE_POSITIVE);

            delete angleCorr;
        }
    }

    // merge result global
    cECorrCUDA.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                                 resDataMaxIndPaddedGlobal.getAddrData(),
                                 XCORR_MERGE_POSITIVE);

    cECorrCUDA.memFree();

    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());
#else
    std::cout << "threadFuncECorrCUDA(): "
              << "this function is not supported without CUDA"
              << std::endl;
#endif
}

#ifndef NDEBUG
void threadFuncNCorr(int threadID, int mpiID)
#else
void threadFuncNCorr(int , int )
#endif
{
    std::stringstream                   strStream;
    boost::posix_time::microseconds     waitingTime(1);
    cData3XCorr<float>                  *angleCorr = NULL;

    // cNCorr object & prepare ref
    cNCorrSingle    cNCorr(mapRefLocal->getNrow(),
                           mapRefLocal->getNcol(),
                           mapRefLocal->getNsec());

    cNCorr.memAlloc();
    cNCorr.prepareRef(mapRefLocal->getAddrData());
    cNCorr.setSizeTpl(mapTplLocal->getNrow(),
                      mapTplLocal->getNcol(),
                      mapTplLocal->getNsec());
    cNCorr.copyTplMskRot3D(mapTplLocal->getAddrData(),
                           mapTplLocal->maskGetAddrData());

    while (streamData || !queueAngle.empty()) {
        if (queueAngle.wait_and_pop(angleCorr, waitingTime)) {
/*#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    rotation"   << std::setw(9) << angleCorr->getID()
                      << ",    threadCPU " << std::setw(6) << threadID
                      << ",    angle1 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(0)
                      << ",    angle2 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(1)
                      << ",    angle3 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(2);
            queueMsge.push(strStream.str());
#endif*/

            // prepare tpl and msk
            cNCorr.tplCreateByRot3D(angleCorr->getElement(0),
                                    angleCorr->getElement(1),
                                    angleCorr->getElement(2),
                                    INTER_LINEAR);
            cNCorr.mskCreateByRot3D(angleCorr->getElement(0),
                                    angleCorr->getElement(1),
                                    angleCorr->getElement(2),
                                    INTER_LINEAR);
            cNCorr.normalizeTpl();
            cNCorr.prepareTplMsk();

            // compute NCC
            cNCorr.computeCorr();

            // merge result
            cNCorr.mergeResult(angleCorr->getID(), XCORR_MERGE_POSITIVE);

            delete angleCorr;
        }
    }

    // merge result global
    cNCorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                             resDataMaxIndPaddedGlobal.getAddrData(),
                             XCORR_MERGE_POSITIVE);

    cNCorr.memFree();
}

#ifndef NDEBUG
void threadFuncNCorrCUDA(int threadID, int mpiID)
#else
#ifdef __GEM_USE_CUDA__
void threadFuncNCorrCUDA(int threadID, int )
#else
void threadFuncNCorrCUDA(int , int )
#endif
#endif
{
#ifdef __GEM_USE_CUDA__
    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());

    std::stringstream                   strStream;
    boost::posix_time::microseconds     waitingTime(1);
    cData3XCorr<float>                  *angleCorr = NULL;

    // cNCorr object & prepare ref
    cuNCorrSingle       cuNCorr(mapRefLocal->getNrow(),
                                mapRefLocal->getNcol(),
                                mapRefLocal->getNsec());

    cuNCorr.memAlloc();
    cuNCorr.prepareRef(mapRefLocal->getAddrData());
    cuNCorr.setSizeTpl(mapTplLocal->getNrow(),
                       mapTplLocal->getNcol(),
                       mapTplLocal->getNsec());
    cuNCorr.setRotTexture(bRotTextureLocal);
    cuNCorr.copyTplMskRot3D(mapTplLocal->getAddrData(),
                            mapTplLocal->maskGetAddrData());

    while (streamData || !queueAngle.empty()) {
        if (queueAngle.wait_and_pop(angleCorr, waitingTime)) {
/*#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    rotation"   << std::setw(9) << angleCorr->getID()
                      << ",    threadGPU " << std::setw(6) << threadID
                      << ",    angle1 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(0)
                      << ",    angle2 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(1)
                      << ",    angle3 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(2);
            queueMsge.push(strStream.str());
#endif*/

            // prepare tpl and msk
            cuNCorr.tplCreateByRot3D(angleCorr->getElement(0),
                                     angleCorr->getElement(1),
                                     angleCorr->getElement(2),
                                     INTER_LINEAR);
            cuNCorr.mskCreateByRot3D(angleCorr->getElement(0),
                                     angleCorr->getElement(1),
                                     angleCorr->getElement(2),
                                     INTER_LINEAR);
            cuNCorr.normalizeTpl();
            cuNCorr.prepareTplMsk();

            // compute NCC
            cuNCorr.computeCorr();

            // merge result
            cuNCorr.mergeResult(angleCorr->getID(), XCORR_MERGE_POSITIVE);

            delete angleCorr;
        }
    }

    // merge result global
    cuNCorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                              resDataMaxIndPaddedGlobal.getAddrData(),
                              XCORR_MERGE_POSITIVE);

    cuNCorr.memFree();

    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());
#else
    std::cout << "threadFuncNCorrCUDA(): "
              << "this function is not supported without CUDA"
              << std::endl;
#endif
}

#ifndef NDEBUG
void threadFuncWCorr(int threadID, int mpiID)
#else
void threadFuncWCorr(int , int )
#endif
{
    std::stringstream                   strStream;
    boost::posix_time::microseconds     waitingTime(1);
    cData3XCorr<float>                  *angleCorr = NULL;

    // cWCorr object & prepare ref
    cWCorrSingle    cWCorr(mapRefLocal->getNrow(),
                           mapRefLocal->getNcol(),
                           mapRefLocal->getNsec());

    cWCorr.memAlloc();
    cWCorr.prepareRef(mapRefLocal->getAddrData());
    cWCorr.setSizeTpl(mapTplLocal->getNrow(),
                      mapTplLocal->getNcol(),
                      mapTplLocal->getNsec());
    cWCorr.copyTplMskWghtRot3D(mapTplLocal->getAddrData(),
                               mapTplLocal->maskGetAddrData(),
                               mapTplLocal->wghtTplGetAddrData(),
                               mapTplLocal->wghtRefGetAddrData());

    while (streamData || !queueAngle.empty()) {
        if (queueAngle.wait_and_pop(angleCorr, waitingTime)) {
/*#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    rotation"   << std::setw(9) << angleCorr->getID()
                      << ",    threadCPU " << std::setw(6) << threadID
                      << ",    angle1 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(0)
                      << ",    angle2 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(1)
                      << ",    angle3 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(2);
            queueMsge.push(strStream.str());
#endif*/

            // prepare tpl and msk
            cWCorr.tplCreateByRot3D(angleCorr->getElement(0),
                                    angleCorr->getElement(1),
                                    angleCorr->getElement(2),
                                    INTER_LINEAR);
            cWCorr.mskCreateByRot3D(angleCorr->getElement(0),
                                    angleCorr->getElement(1),
                                    angleCorr->getElement(2),
                                    INTER_LINEAR);
            cWCorr.wghtTplCreateByRot3D(angleCorr->getElement(0),
                                        angleCorr->getElement(1),
                                        angleCorr->getElement(2),
                                        INTER_LINEAR);
            cWCorr.wghtRefCreateByRot3D(angleCorr->getElement(0),
                                        angleCorr->getElement(1),
                                        angleCorr->getElement(2),
                                        INTER_LINEAR);
            cWCorr.normalizeTpl();
            cWCorr.prepareTplRefWght();

            // compute WCC
            cWCorr.computeCorr();

            // merge result
            cWCorr.mergeResult(angleCorr->getID(), XCORR_MERGE_POSITIVE);

            delete angleCorr;
        }
    }

    // merge result global
    cWCorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                             resDataMaxIndPaddedGlobal.getAddrData(),
                             XCORR_MERGE_POSITIVE);

    cWCorr.memFree();
}

#ifndef NDEBUG
void threadFuncWCorrCUDA(int threadID, int mpiID)
#else
#ifdef __GEM_USE_CUDA__
void threadFuncWCorrCUDA(int threadID, int )
#else
void threadFuncWCorrCUDA(int , int )
#endif
#endif
{
#ifdef __GEM_USE_CUDA__
    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());

    std::stringstream                   strStream;
    boost::posix_time::microseconds     waitingTime(1);
    cData3XCorr<float>                  *angleCorr = NULL;

    // cWCorr object & prepare ref
    cuWCorrSingle       cuWCorr(mapRefLocal->getNrow(),
                                mapRefLocal->getNcol(),
                                mapRefLocal->getNsec());

    cuWCorr.memAlloc();
    cuWCorr.prepareRef(mapRefLocal->getAddrData());
    cuWCorr.setSizeTpl(mapTplLocal->getNrow(),
                       mapTplLocal->getNcol(),
                       mapTplLocal->getNsec());
    cuWCorr.setRotTexture(bRotTextureLocal);
    cuWCorr.copyTplMskWghtRot3D(mapTplLocal->getAddrData(),
                                mapTplLocal->maskGetAddrData(),
                                mapTplLocal->wghtTplGetAddrData(),
                                mapTplLocal->wghtRefGetAddrData());

    while (streamData || !queueAngle.empty()) {
        if (queueAngle.wait_and_pop(angleCorr, waitingTime)) {
/*#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    rotation"   << std::setw(9) << angleCorr->getID()
                      << ",    threadGPU " << std::setw(6) << threadID
                      << ",    angle1 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(0)
                      << ",    angle2 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(1)
                      << ",    angle3 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(2);
            queueMsge.push(strStream.str());
#endif*/

            // prepare tpl and msk
            cuWCorr.tplCreateByRot3D(angleCorr->getElement(0),
                                     angleCorr->getElement(1),
                                     angleCorr->getElement(2),
                                     INTER_LINEAR);
            cuWCorr.mskCreateByRot3D(angleCorr->getElement(0),
                                     angleCorr->getElement(1),
                                     angleCorr->getElement(2),
                                     INTER_LINEAR);
            cuWCorr.wghtTplCreateByRot3D(angleCorr->getElement(0),
                                         angleCorr->getElement(1),
                                         angleCorr->getElement(2),
                                         INTER_LINEAR);
            cuWCorr.wghtRefCreateByRot3D(angleCorr->getElement(0),
                                         angleCorr->getElement(1),
                                         angleCorr->getElement(2),
                                         INTER_LINEAR);
            cuWCorr.normalizeTpl();
            cuWCorr.prepareTplRefWght();

            // compute WCC
            cuWCorr.computeCorr();

            // merge result
            cuWCorr.mergeResult(angleCorr->getID(), XCORR_MERGE_POSITIVE);

            delete angleCorr;
        }
    }

    // merge result global
    cuWCorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                              resDataMaxIndPaddedGlobal.getAddrData(),
                              XCORR_MERGE_POSITIVE);

    cuWCorr.memFree();

    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());
#else
    std::cout << "threadFuncWCorrCUDA(): "
              << "this function is not supported without CUDA"
              << std::endl;
#endif
}

#ifndef NDEBUG
void threadFuncCCorr(int threadID, int mpiID)
#else
void threadFuncCCorr(int , int )
#endif
{
    std::stringstream                   strStream;
    boost::posix_time::microseconds     waitingTime(1);
    cData3XCorr<float>                  *angleCorr = NULL;

    // cCCorr object & prepare ref
    cCCorrSingle       cCCorr(mapRefLocal->getNrow(),
                              mapRefLocal->getNcol(),
                              mapRefLocal->getNsec());

    cCCorr.memAlloc();

    cCCorr.getRefCorrDensity()->prepareRef(mapRefLocal->getAddrData());
    cCCorr.getRefCorrSurface()->prepareRef(mapRefSurfaceLocal->getAddrData());
    cCCorr.getRefCorrPenalty()->prepareRef(mapRefOuterLocal->getAddrData());

    cCCorr.setSizeTpl(mapTplLocal->getNrow(),
                      mapTplLocal->getNcol(),
                      mapTplLocal->getNsec());
    cCCorr.setWeight(wDensityLocal, wSurfaceLocal, wPenaltyLocal);

    cCCorr.getRefCorrDensity()->copyTplMskWghtRot3D(mapTplLocal->getAddrData(),
                                                    mapTplLocal->maskGetAddrData(),
                                                    mapTplLocal->wghtTplGetAddrData(),
                                                    mapTplLocal->wghtRefGetAddrData());
    cCCorr.getRefCorrSurface()->copyTplRot3D(mapTplSurfaceLocal->getAddrData());
    cCCorr.getRefCorrPenalty()->copyTplRot3D(mapTplInnerLocal->getAddrData());

    while (streamData || !queueAngle.empty()) {
        if (queueAngle.wait_and_pop(angleCorr, waitingTime)) {

/*#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    rotation"   << std::setw(9) << angleCorr->getID()
                      << ",    threadGPU " << std::setw(6) << threadID
                      << ",    angle1 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(0)
                      << ",    angle2 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(1)
                      << ",    angle3 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(2);
            queueMsge.push(strStream.str());
#endif*/

            // prepare tpl and msk
            cCCorr.rotateTpl(angleCorr->getElement(0),
                             angleCorr->getElement(1),
                             angleCorr->getElement(2),
                             INTER_LINEAR);
            cCCorr.normalizeTpl();
            cCCorr.prepareTpl();

            // compute CCC
            cCCorr.computeCorr();

            // merge result
            cCCorr.mergeResult(angleCorr->getID(), XCORR_MERGE_NEGATIVE);   // negative for surface fitting

            delete angleCorr;
        }
    }

    // merge result global
    cCCorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                              resDataMaxIndPaddedGlobal.getAddrData(),
                              XCORR_MERGE_POSITIVE);

    cCCorr.memFree();
}

#ifndef NDEBUG
void threadFuncCCorrCUDA(int threadID, int mpiID)
#else
#ifdef __GEM_USE_CUDA__
void threadFuncCCorrCUDA(int threadID, int )
#else
void threadFuncCCorrCUDA(int , int )
#endif
#endif
{
#ifdef __GEM_USE_CUDA__
    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());

    std::stringstream                   strStream;
    boost::posix_time::microseconds     waitingTime(1);
    cData3XCorr<float>                  *angleCorr = NULL;

    // cCCorr object & prepare ref
    cuCCorrSingle       cuCCorr(mapRefLocal->getNrow(),
                                mapRefLocal->getNcol(),
                                mapRefLocal->getNsec());

    cuCCorr.memAlloc();

    cuCCorr.getRefCorrDensity()->prepareRef(mapRefLocal->getAddrData());
    cuCCorr.getRefCorrSurface()->prepareRef(mapRefSurfaceLocal->getAddrData());
    cuCCorr.getRefCorrPenalty()->prepareRef(mapRefOuterLocal->getAddrData());

    cuCCorr.setSizeTpl(mapTplLocal->getNrow(),
                       mapTplLocal->getNcol(),
                       mapTplLocal->getNsec());
    cuCCorr.setRotTexture(bRotTextureLocal);
    cuCCorr.setWeight(wDensityLocal, wSurfaceLocal, wPenaltyLocal);

    cuCCorr.getRefCorrDensity()->copyTplMskWghtRot3D(mapTplLocal->getAddrData(),
                                                     mapTplLocal->maskGetAddrData(),
                                                     mapTplLocal->wghtTplGetAddrData(),
                                                     mapTplLocal->wghtRefGetAddrData());
    cuCCorr.getRefCorrSurface()->copyTplRot3D(mapTplSurfaceLocal->getAddrData());
    cuCCorr.getRefCorrPenalty()->copyTplRot3D(mapTplInnerLocal->getAddrData());

    while (streamData || !queueAngle.empty()) {
        if (queueAngle.wait_and_pop(angleCorr, waitingTime)) {

/*#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    rotation"   << std::setw(9) << angleCorr->getID()
                      << ",    threadGPU " << std::setw(6) << threadID
                      << ",    angle1 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(0)
                      << ",    angle2 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(1)
                      << ",    angle3 = "  << std::setprecision(6)
                      << std::setw(8) << std::fixed << angleCorr->getElement(2);
            queueMsge.push(strStream.str());
#endif*/

            // prepare tpl and msk
            cuCCorr.rotateTpl(angleCorr->getElement(0),
                              angleCorr->getElement(1),
                              angleCorr->getElement(2),
                              INTER_LINEAR);
            cuCCorr.normalizeTpl();
            cuCCorr.prepareTpl();

            // compute CCC
            cuCCorr.computeCorr();

            // merge result
            cuCCorr.mergeResult(angleCorr->getID(), XCORR_MERGE_NEGATIVE);  // negative for surface fitting

            delete angleCorr;
        }
    }

    // merge result global
    cuCCorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                              resDataMaxIndPaddedGlobal.getAddrData(),
                              XCORR_MERGE_POSITIVE);

    cuCCorr.memFree();

    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());
#else
    std::cout << "threadFuncCCorrCUDA(): "
              << "this function is not supported without CUDA"
              << std::endl;
#endif
}

void threadFuncData(int mpiID, int mpiNumProcs,
                    int nNode, int nRotZ)
{
    std::stringstream           strStream;
    size_t                      nRotation = (size_t) (nNode*nRotZ);
    cVector3<float>             angle;
    cData3XCorr<float>          *angleCorr = NULL;
    int                         percentRead = -1;

    //for (size_t iRot = mpiID; iRot < 1; iRot += mpiNumProcs) {
    for (size_t iRot = mpiID; iRot < nRotation; iRot += mpiNumProcs) {
        if (mpiID == 0) {
            if ((int) (100*(iRot+1)/nRotation) > percentRead) {
                percentRead = (int) (100*(iRot+1)/nRotation);
            }
        }

#ifndef NDEBUG
        // message
        strStream.str(std::string());
        strStream << "Process" << std::setw(5) << mpiID
                  << ",    rotation" << std::setw(9) << iRot+1
                  << ",    iNode = " << std::setw(5) << iRot / nRotZ
                  << ",    iRotZ = " << iRot - (iRot / nRotZ) * nRotZ;
        queueMsge.push(strStream.str());
#else
        if (mpiID == 0) {
            std::cout << "\r     " << std::setw(3) << percentRead << "%     " << std::flush;
        }
#endif
        if (mpiID == 0) {
            std::cout << "\r     " << std::setw(3) << percentRead << "%     " << std::flush;
        }

        pft_angle(iRot, nNode, nRotZ, angle);

        // new angle
        angleCorr = new cData3XCorr<float>(cVector3<size_t>(3,1,1));
        angleCorr->setElement(0, angle[0]);
        angleCorr->setElement(1, angle[1]);
        angleCorr->setElement(2, angle[2]);
        angleCorr->setID(iRot+1);

        // push to queue
        queueAngle.wait_and_push(angleCorr);
    }

    streamData = false;
}

void threadFuncMsge(void)
{
    std::string                         strMsge;
    boost::posix_time::milliseconds     sleepingTime(1);

    while (streamMsge || !queueMsge.empty()) {
        if (queueMsge.try_pop(strMsge)) {
            std::cout << strMsge << std::endl;
        }
        else {
            boost::this_thread::sleep(sleepingTime);
        }
    }
}

void cFitter::computeCorr(void)
{
    const std::string   funcName("void cFitter::computeCorr(void)");

    bRotTextureLocal = bRotTexture;

    wDensityLocal = (float) wDensity;
    wSurfaceLocal = (float) wSurface;
    wPenaltyLocal = (float) wPenalty;

    mapRefLocal      = &mapRef;
    mapTplLocal      = &mapTpl;

    mapRefSurfaceLocal = &mapRefSurface;
    mapTplSurfaceLocal = &mapTplSurface;
    mapRefOuterLocal   = &mapRefOuter;
    mapTplInnerLocal   = &mapTplInner;

    // RESULT: allocate
    resDataAbsMaxPaddedGlobal.memReAllocZero(mapRefLocal->getSize());
    resDataMaxIndPaddedGlobal.memReAlloc    (mapRefLocal->getSize());

#ifdef __GEM_USE_MPI__
        if (mpiNumProcs > 1) { MPI_Barrier(MPI_COMM_WORLD); }
#endif

    // THREAD: init
    streamData = true;
    streamMsge = true;

    boost::thread   threadData(threadFuncData, mpiID, mpiNumProcs, nNode, nRotZ);
    boost::thread   threadMsge(threadFuncMsge);
    boost::thread   **threadCompute;

    // THREAD: allocate
    threadCompute = new boost::thread *[nGPU+nCPU];
    for (int i = 0; i < nGPU+nCPU; i++) {
        if (i < nGPU) {
            switch (mCorr) {
                case CC_MODE_XCORR: threadCompute[i] = new boost::thread(threadFuncXCorrCUDA, i+1, mpiID); break;
                case CC_MODE_ECORR: threadCompute[i] = new boost::thread(threadFuncECorrCUDA, i+1, mpiID); break;
                case CC_MODE_NCORR: threadCompute[i] = new boost::thread(threadFuncNCorrCUDA, i+1, mpiID); break;
                case CC_MODE_WCORR: threadCompute[i] = new boost::thread(threadFuncWCorrCUDA, i+1, mpiID); break;
                case CC_MODE_SCORR: threadCompute[i] = new boost::thread(threadFuncXCorrCUDA, i+1, mpiID); break;
                case CC_MODE_CCORR: threadCompute[i] = new boost::thread(threadFuncCCorrCUDA, i+1, mpiID); break;
                default:            ERROR(funcName, "unsupported correlation mode");
            }
        }
        else {
            switch (mCorr) {
                case CC_MODE_XCORR: threadCompute[i] = new boost::thread(threadFuncXCorr, i+1, mpiID); break;
                case CC_MODE_ECORR: threadCompute[i] = new boost::thread(threadFuncECorr, i+1, mpiID); break;
                case CC_MODE_NCORR: threadCompute[i] = new boost::thread(threadFuncNCorr, i+1, mpiID); break;
                case CC_MODE_WCORR: threadCompute[i] = new boost::thread(threadFuncWCorr, i+1, mpiID); break;
                case CC_MODE_SCORR: threadCompute[i] = new boost::thread(threadFuncXCorr, i+1, mpiID); break;
                case CC_MODE_CCORR: threadCompute[i] = new boost::thread(threadFuncCCorr, i+1, mpiID); break;
                default:            ERROR(funcName, "unsupported correlation mode");
            }
        }
    }

    // THREAD: join
    threadData.join();
    for (int i = 0; i < nGPU+nCPU; i++) { threadCompute[i]->join(); }
    streamMsge = false;
    threadMsge.join();

    // THREAD: deallocate
    for (int i = 0; i < nGPU+nCPU; i++) { delete threadCompute[i]; }
    delete [] threadCompute;

#ifdef __GEM_USE_MPI__
    if (mpiNumProcs > 1) { MPI_Barrier(MPI_COMM_WORLD); }
#endif

    // RESULT: reduction
    if (mpiNumProcs > 1) {
#ifdef __GEM_USE_MPI__

#ifdef __GEM_USE_OPENMP__
        omp_set_num_threads(omp_get_num_procs());
#endif  // __GEM_USE_OPENMP__

        mpiReduce_pickerV2(resDataAbsMaxPaddedGlobal.getAddrData(),
                           resDataMaxIndPaddedGlobal.getAddrData(),
                           mapRefLocal->getNelement(),
                           XCORR_MERGE_POSITIVE);

#ifdef __GEM_USE_OPENMP__
        omp_set_num_threads(1);
#endif  // __GEM_USE_OPENMP__

#endif  // __GEM_USE_MPI__
    }

    // RESULT: save and deallocate
    switch (mCorr) {
        case CC_MODE_XCORR:
            break;
        case CC_MODE_ECORR:     // remove values not in the range [-1 1]
            resDataAbsMaxPaddedGlobal.opLimitCorrValue(-1,1);
            break;
        case CC_MODE_NCORR:     // remove values not in the range [-1 1]
            resDataAbsMaxPaddedGlobal.opLimitCorrValue(-1,1);
            break;
        case CC_MODE_WCORR:     // remove values not in the range [-1 1]
            resDataAbsMaxPaddedGlobal.opLimitCorrValue(-1,1);
            break;
        case CC_MODE_SCORR:
            break;
        case CC_MODE_CCORR:
            break;
        default:
            ERROR(funcName, "unsupported correlation mode");
    }

    cVector3<ptrdiff_t>     offset((mapTplLocal->getNrow()-1)/2,
                                   (mapTplLocal->getNcol()-1)/2,
                                   (mapTplLocal->getNsec()-1)/2);

    resDataAbsMaxPaddedGlobal.opCircShift(offset);
    resDataMaxIndPaddedGlobal.opCircShift(offset);

    mapCor.opTypeCast(resDataAbsMaxPaddedGlobal);
    mapInd.opTypeCast(resDataMaxIndPaddedGlobal);

    mapCor.opFFTSizeRestore(sizeRefOriginal);
    mapInd.opFFTSizeRestore(sizeRefOriginal);

    resDataAbsMaxPaddedGlobal.memFree();
    resDataMaxIndPaddedGlobal.memFree();
}

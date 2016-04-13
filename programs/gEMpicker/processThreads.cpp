/***********************************************************************
 *  File:       processThreads.cpp
 *
 *  Purpose:    Implementation of threading functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "gEMpicker.hpp"

static cMapHeader                               mapRefHeader;
static cData3XCorr<float>                       *mapRef = NULL;
static cQueueConcurrent<cData3XCorr<float>* >   queueData(30);
static cQueueConcurrent<std::string>            queueMsge;
static bool                                     streamData, streamMsge;
static bool                                     bMskSingle;

static cArray3<float>       resDataAbsMaxPaddedGlobal;
static cArray3<size_t>      resDataMaxIndPaddedGlobal;
static cVector3<size_t>     tplSize;

size_t                      iRefStep, iTplStep;

#ifndef NDEBUG
void threadFuncNCorr(int threadID, int mpiID)
#else
void threadFuncNCorr(int , int )
#endif
{
    std::stringstream                    strStream;
    boost::posix_time::microseconds      waitingTime(1);
    cData3XCorr<float>                   *mapFromQueue = NULL;

    // cNCorr object & prepare ref
    cNCorrSingle    cNCorr(mapRef->getNrow(),
                           mapRef->getNcol());

    cNCorr.memAlloc();
    cNCorr.prepareRef(mapRef->getAddrData());
    if (bMskSingle) cNCorr.setMskSingle();

    while (streamData || !queueData.empty()) {
        if (queueData.wait_and_pop(mapFromQueue, waitingTime)) {
#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    thread " << std::setw(7) << threadID
                      << ",    CPU" << std::setw(7)
                      << mapFromQueue->getID();
            queueMsge.push(strStream.str());
#endif

            // prepare tpl and msk
            cNCorr.setSizeTpl(mapFromQueue->getNrow(),
                              mapFromQueue->getNrow());
            cNCorr.copyTplMsk(mapFromQueue->getAddrData(),
                              mapFromQueue->maskGetAddrData());
            cNCorr.normalizeTpl(mapFromQueue->maskGetArea());
            cNCorr.prepareTplMsk();

            // compute NCC
            cNCorr.computeCorr();

            // merge result
            cNCorr.mergeResult(mapFromQueue->getID(), XCORR_MERGE_NEGATIVE);

            delete mapFromQueue;
        }
    }

    // merge result global
    cNCorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                             resDataMaxIndPaddedGlobal.getAddrData(),
                             XCORR_MERGE_NEGATIVE);

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

    std::stringstream                    strStream;
    boost::posix_time::microseconds      waitingTime(1);
    cData3XCorr<float>                   *mapFromQueue = NULL;

    // cNCorr object & prepare ref
    cuNCorrSingle       cuNCorr(mapRef->getNrow(),
                                mapRef->getNcol());

    cuNCorr.memAlloc();
    cuNCorr.prepareRef(mapRef->getAddrData());
    if (bMskSingle) cuNCorr.setMskSingle();

    while (streamData || !queueData.empty()) {
        if (queueData.wait_and_pop(mapFromQueue, waitingTime)) {
#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    thread " << std::setw(7) << threadID
                      << ",    GPU" << std::setw(7)
                      << mapFromQueue->getID();
            queueMsge.push(strStream.str());
#endif

            // prepare tpl and msk
            cuNCorr.setSizeTpl(mapFromQueue->getNrow(),
                               mapFromQueue->getNcol());
#ifdef __GEM_USE_CUDA_PINNEDMEM__
            cuNCorr.copyTplMsk(mapFromQueue->getRefDataPinned(),
                               mapFromQueue->getRefMaskPinned());
#else
            cuNCorr.copyTplMsk(mapFromQueue->getAddrData(),
                               mapFromQueue->maskGetAddrData());
#endif
            cuNCorr.normalizeTpl(mapFromQueue->maskGetArea());
            cuNCorr.prepareTplMsk();

            // compute NCC
            cuNCorr.computeCorr();

            // merge result
            cuNCorr.mergeResult(mapFromQueue->getID(), XCORR_MERGE_NEGATIVE);

            delete mapFromQueue;
        }
    }

    // merge result global
    cuNCorr.mergeResultGlobal(resDataAbsMaxPaddedGlobal.getAddrData(),
                              resDataMaxIndPaddedGlobal.getAddrData(),
                              XCORR_MERGE_NEGATIVE);

    cuNCorr.memFree();

    CUDA_SAFE_CALL(cudaSetDevice(threadID-1));
    CUDA_SAFE_CALL(cudaDeviceReset());
#else
    std::cout << "threadFuncNCorrCUDA(): "
              << "this function is not supported without CUDA"
              << std::endl;
#endif
}

void threadFuncData(int mpiID,
                    float angle2D, size_t numRot2D,
                    std::vector<std::string>& fileListTpl,
                    std::vector<std::string>& fileListMsk)
{
    std::stringstream           strStream;
    cData3XCorr<float>          tplFromDisk, mskFromDisk, *mapRotated;
    int                         percentRead = -1;
    size_t                      indx;
    cVector3<float>             angle = 0.0f;

    // process each template
    for (size_t iTpl = 0; iTpl < fileListTpl.size(); iTpl += iTplStep) {
        // read a template form disk
        cMRC<float>().read(fileListTpl[iTpl], tplFromDisk);
        if (iTpl == (size_t) mpiID) tplSize  = tplFromDisk.getSize();

        // read a mask form disk
        if (bMskSingle) {
            if (mskFromDisk.getSize() == (size_t) 0)  cTIFF<float>().read(fileListMsk[0], mskFromDisk);
        }
        else {
            cTIFF<float>().read(fileListMsk[iTpl], mskFromDisk);
        }

        // conditions on dimension and size
        require(tplFromDisk.getSize() <= mapRef->getSize(), "search data is too-large");

        require(tplFromDisk.getDimension() == 2,  "search data dimension != 2");
        require(tplFromDisk.getSize() == tplSize, "search data have different size");

        require(mskFromDisk.getDimension() == 2,  "mask data dimension != 2");
        require(mskFromDisk.getSize() == tplSize, "mask data have different size");

        // scan all rotations
        for (size_t iRot = 0; iRot < numRot2D; iRot++) {
            indx = iTpl*numRot2D + iRot;

            if (mpiID == 0) {
                if ((int) (100*(indx+1)/fileListTpl.size()) > percentRead) {
                    percentRead = (int) (100*(indx+1)/(numRot2D*fileListTpl.size()));
                }
            }

#ifndef NDEBUG
            // message
            strStream.str(std::string());
            strStream << "Process" << std::setw(5) << mpiID
                      << ",    pushing" << std::setw(7) << indx+1
                      << ",    " << fileListTpl[iTpl]
                      << ",    angle" << std::setw(7) << (float) iRot*angle2D
                      << ",    ";
            if (bMskSingle) strStream << fileListMsk[0];
            else            strStream << fileListMsk[iTpl];
            queueMsge.push(strStream.str());
#else
            if (mpiID == 0) {
                std::cout << "\r     " << std::setw(3) << percentRead << "%     " << std::flush;
            }
#endif

            angle[0] = deg2rad((float) iRot * angle2D);

            // new rotation
            mapRotated = new cData3XCorr<float>;

            cTransSimilarity<float>().rotate(tplFromDisk, *mapRotated, angle);

            if (bMskSingle) { mapRotated->maskCreateByCopy(mskFromDisk);          }
            else            { mapRotated->maskCreateByRotation(mskFromDisk, angle); }

            // set indx
            mapRotated->setID(indx+1);

            // normalize the template
            mapRotated->opNormalizeTplNCC(true);

            // prepare pinned memory
#ifdef __GEM_USE_CUDA_PINNEDMEM__
            mapRotated->preparePinnedMem();
#endif

            // push to queue
            queueData.wait_and_push(mapRotated);
        }
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

void threadFuncInit(sPickerParams& params)
{
    std::stringstream   strStream;
    std::ofstream       outfile;
    cTime               cTimer;
    double              timeComputeTotal = 0;
    cVector3<size_t>    sizeRefOri;

    bMskSingle = params.bMskSingle;

    mapRef = new cData3XCorr<float>;

    // control the distribution mode
    if (params.mpiDistTgt) {
        iRefStep = params.mpiNumProcs;
        iTplStep = 1;
    }
    else {
        iRefStep = 1;
        iTplStep = params.mpiNumProcs;
    }

    // process each reference
    //for (size_t iRef = 0; iRef < 1; iRef += iRefStep) {
    for (size_t iRef = 0; iRef < params.fileListRef.size(); iRef += iRefStep) {
        if (params.mpiID == 0) {
            std::cout << "Processing micrograph "
                      << iRef+1 << "/" << params.fileListRef.size() << ": "
                      << params.fileListRef[iRef] << std::endl;
        }

        // TIME_COMPUTE
        cTimer.tic();

        // REFERENCE: allocate
        cMRC<float>().read(params.fileListRef[iRef], *mapRef, &mapRefHeader);
        require(mapRef->getDimension() == 2,
                "micrograph dimension != 2");
        mapRef->opFFTSizeChange(sizeRefOri, cVector3<size_t>(0,0,0));

#ifdef __GEM_LIMIT_FUNC__
        require(mapRef->getNrow() <= 4096, "micrograph data is too-large");
        require(mapRef->getNcol() <= 4096, "micrograph data is too-large");
#endif

        // RESULT: allocate
        resDataAbsMaxPaddedGlobal.memReAllocZero(mapRef->getSize());
        resDataMaxIndPaddedGlobal.memReAlloc    (mapRef->getSize());

        // THREAD: init
        streamData = true;
        streamMsge = true;

        boost::thread   threadData(threadFuncData,
                                   params.mpiID,
                                   params.angle2D,
                                   params.numRot2D,
                                   params.fileListTpl,
                                   params.fileListMsk);
        boost::thread   threadMsge(threadFuncMsge);
        boost::thread   **threadCompute;

        // THREAD: allocate
        threadCompute = new boost::thread *[params.nGPU+params.nCPU];
        for (int i = 0; i < params.nGPU+params.nCPU; i++) {
            if (i < params.nGPU) {
                threadCompute[i] = new boost::thread(
                                            threadFuncNCorrCUDA,
                                            i+1, params.mpiID);
            }
            else {
                threadCompute[i] = new boost::thread(
                                            threadFuncNCorr,
                                            i+1, params.mpiID);
            }
        }

        // THREAD: join
        threadData.join();
        for (int i = 0; i < params.nGPU+params.nCPU; i++) {
            threadCompute[i]->join();
        }
        streamMsge = false;
        threadMsge.join();

        // THREAD: deallocate
        for (int i = 0; i < params.nGPU+params.nCPU; i++) {
            delete threadCompute[i];
        }
        delete [] threadCompute;

#ifdef __GEM_USE_MPI__
        if (params.mpiNumProcs > 1) {
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif

        // TIME_COMPUTE
        cTimer.toc("noshow");
        timeComputeTotal += cTimer.getElapsedTime();

        if (params.mpiID == 0) {
            strStream << "Correlation time:" << std::fixed
                      << std::setw(5) << params.mpiNumProcs << " nProc, "
                      << std::setw(5) << params.nGPU << " nGPU, "
                      << std::setw(5) << params.nCPU << " nCPU, "
                      << "   " << cTimer.getElapsedTime()
                      << std::endl;

#ifdef __GEM_SAVE_TIMING__
            outfile.open("time_gempicker_compute.txt", std::ofstream::app |
                                                       std::ofstream::out);
            outfile << strStream.str();
            outfile.close();
#endif  // __GEM_SAVE_TIMING__

            std::cout << strStream.str();
            strStream.str("");
        }

        // RESULT: reduction
        if (params.mpiNumProcs > 1 && !params.mpiDistTgt) {
            // TIME_REDUCTION
            cTimer.tic();

#ifdef __GEM_USE_MPI__

#ifdef __GEM_USE_OPENMP__
            omp_set_num_threads(omp_get_num_procs());
#endif  // __GEM_USE_OPENMP__

            mpiReduce_pickerV2(resDataAbsMaxPaddedGlobal.getAddrData(),
                               resDataMaxIndPaddedGlobal.getAddrData(),
                               mapRef->getNelement(),
                               XCORR_MERGE_NEGATIVE);
            /* This code is for comparison purpose only
            mpiReduce_pickerV1(resDataAbsMaxPaddedGlobal.getAddrData(),
                               resDataMaxIndPaddedGlobal.getAddrData(),
                               mapRef->getNelement(),
                               XCORR_MERGE_NEGATIVE);
            */

#ifdef __GEM_USE_OPENMP__
            omp_set_num_threads(1);
#endif  // __GEM_USE_OPENMP__

#endif  // __GEM_USE_MPI__

            // TIME_REDUCTION
            cTimer.toc("noshow");

            if (params.mpiID == 0) {
                strStream << "Reduction time:  " << std::fixed
                          << std::setw(5)    << params.mpiNumProcs << " nProc, "
                          << "   " << cTimer.getElapsedTime()
                          << std::endl;

#ifdef __GEM_SAVE_TIMING__
                outfile.open("time_gempicker_reduction.txt", std::ofstream::app |
                                                             std::ofstream::out);
                outfile << strStream.str();
                outfile.close();
#endif  // __GEM_SAVE_TIMING__

                std::cout << strStream.str();
                strStream.str("");
            }
        }

        // RESULT: save and deallocate
        std::string             fileNameRoot = cFileSystem(params.fileListRef[iRef]).getFileRoot();
        std::string             fileNameSave;
        cVector3<ptrdiff_t>     offset((tplSize[0]-1)/2,
                                       (tplSize[1]-1)/2,
                                       (tplSize[2]-1)/2);

        resDataAbsMaxPaddedGlobal.opCircShift(offset);
        resDataMaxIndPaddedGlobal.opCircShift(offset);

        if ((!params.mpiDistTgt && params.mpiID == 0) ||
              params.mpiDistTgt) {
            // TIME_SAVING
            cTimer.tic();
            mapRefHeader.data_mode = 2; // for saving corr results in float32

            cData3XCorr<float>   corrToDisk;
            corrToDisk.opTypeCast(resDataAbsMaxPaddedGlobal);
            corrToDisk.opFFTSizeRestore(sizeRefOri);
            fileNameSave = params.dirResStr + "xcorr/" + fileNameRoot + "__corr.mrc";
            cMRC<float>().write(corrToDisk, fileNameSave, &mapRefHeader);

            fileNameSave = params.dirResStr + "xcorr/" + fileNameRoot + "__corr.tif";
#ifdef __GEM_USE_OPENCV__
            std::string     textStr = num2str(corrToDisk.getMin(), 4) + " : " +
                                      num2str(corrToDisk.getMax(), 4);
            cOpenCV<float>().putText(corrToDisk, fileNameSave, textStr, COLOR_BLUE,
                                     3.0f*(float)corrToDisk.getNrow()/2048.0f, NORM8BIT_TRUE);
#else
            cTIFF<float>().write(corrToDisk, fileNameSave, NORM8BIT_TRUE);
#endif

            cData3XCorr<float>   indxToDisk;
            indxToDisk.opTypeCast(resDataMaxIndPaddedGlobal);
            indxToDisk.opFFTSizeRestore(sizeRefOri);
            fileNameSave = params.dirResStr + "xcorr/" + fileNameRoot + "__indx.mrc";
            cMRC<float>().write(indxToDisk, fileNameSave, &mapRefHeader);

            // TIME_SAVING
            cTimer.toc("noshow");

            std::cout << std::endl;
        }
    }

    resDataAbsMaxPaddedGlobal.memFree();
    resDataMaxIndPaddedGlobal.memFree();

    delete mapRef;

    if (params.mpiID == 0) {
        std::cout << "Total correlation time: " << timeComputeTotal << std::endl;
    }
}

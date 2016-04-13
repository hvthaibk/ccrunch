/***********************************************************************
 *  File:       cAlignerCurtain.cpp
 *
 *  Purpose:    Implementation of cAligner class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2014 Thai V. Hoang, IGBMC
 **********************************************************************/

#include "cAligner.hpp"

// hoang-Mint
// reading time: elapsed time: 31.173
// reading + writing: elapsed time: 135.72

// microgpu2
// reading time: elapsed time: 3.11995
// reading + writing: elapsed time: 111.897

// parameters used in cResultAligner::checkResultRigid()
//  - radius
//  - sigmaNotch
//  - sigmaLowpass

void cAligner::removeCurtain(void)
{
    // output directory
    dirDataTmp1 = dirDataIn;
    dirDataTmp2 = dirDataOut + "/filter";

    cFileSystem     objFileSystem(dirDataTmp2);

    if (objFileSystem.isExist()) {
        objFileSystem.dirDelete();
    }
    objFileSystem.dirCreate();

    // variables
    cData3<float>               image, kernel;
    cFilterFrequency<float>     objFilter;
    cTIFF<float>                objTIFF;
    cTime                       timer;

#ifdef __GEM_USE_CUDA__
    cuData3<float>              imageCUDA;
    cuFilterFrequency<float>    objFilterCUDA;
#endif

    // define the composite kernel
    objTIFF.read(dirDataTmp1 + "/" + fileList[0], image);

    objFilter.kernelFreq2NotchLineHor(image.getSize(), 10, 2);
    kernel  = objFilter.getKernel();

    objFilter.kernelFreq2LowPass(image.getSize(), (float) image.getNrow()/10.f);
    kernel *= objFilter.getKernel();

    objFilter.setKernel(kernel);

    // write kernel to disk
    objTIFF.write(objFilter.getKernelFull2(image.getSize()),
                  dirDataOut + "/kernelFreq.tif", NORM8BIT_TRUE);

#ifdef __GEM_USE_CUDA__
    objFilterCUDA.copyKernelToGPU(objFilter);
#endif

    // frequency filtering
    mpiMasterPrint(mpiID, "Removing curtains...\n");
    timer.tic();
    for (size_t i = 0; i < fileList.size(); i++) {

#ifndef NDEBUG
        std::cout << fileList[i] << std::endl;
#else
        std::cout << "\r     " << i+1 << "/" << num2str(fileList.size()) << std::flush;
#endif

        objTIFF.read(dirDataTmp1 + "/" + fileList[i], image);

#ifdef __GEM_USE_CUDA__
        imageCUDA = image;
        objFilterCUDA.filterFreq2(imageCUDA);
        image = imageCUDA;
#else
        objFilter.filterFreq2(image);
#endif

        objTIFF.write(image, dirDataTmp2 + "/" + fileList[i]);
    }
    timer.toc("     filtering time: ");
    std::cout << std::endl;
}

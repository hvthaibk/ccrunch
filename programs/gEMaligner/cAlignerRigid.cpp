/***********************************************************************
 *  File:       cAlignerTranslation.cpp
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

// parameters used in cResultAligner::checkResultRigid()
//  - thrOffset
//  - thrCorrVal

void cAligner::alignRigid(void)
{
    // output directory
    dirDataTmp1 = dirDataOut + "/filter";
    dirDataTmp2 = dirDataOut + "/rigid";

    cFileSystem     objFileSystem(dirDataTmp2);

    if (objFileSystem.isExist()) {
        objFileSystem.dirDelete();
    }
    objFileSystem.dirCreate();

    // variables
    cData3<float>                   image;
    cSize3                          imSize;
    cStackRegTrans<float>           objStackRegTrans;
    std::vector< cVector3<float> >  corrRes;
    cTime                           timer;

    // ---------------------
    // Define ROIs
    // ---------------------
    cTIFF<float>().read(dirDataTmp1 + "/" + fileList[0], image);
    imSize = image.getSize();

    objStackRegTrans.setROIs(cSize3(imSize[0]/2,imSize[1]/2,1),
                             cSize3(imSize[0]/4,imSize[1]/4,0),
                             imSize,
                             cSize3(0,0,0));

    // ---------------------
    // Compute offsets
    // ---------------------
    startNthreadsFFTW();

    mpiMasterPrint(mpiID, "Computing displacements between images in stack...\n");
    timer.tic();
    objStackRegTrans.computeOffsetsUsingCorr(dirDataTmp1, fileList);
    timer.toc("     computation time: ");
    std::cout << std::endl;

    stopNthreadsFFTW();

    // ---------------------
    // Check offsets
    // ---------------------
    mpiMasterPrint(mpiID, "Checking correlation results...\n");
    corrRes = objStackRegTrans.getResult();
    cResultAligner<float>().checkResultRigid(corrRes, imSize, 0.1f, 0.9f);
    std::cout << std::endl;

    // ---------------------
    // Align sections
    // ---------------------
    startOpenMP();

    mpiMasterPrint(mpiID, "Aligning images...\n");
    timer.tic();
    objStackRegTrans.alignSections(dirDataTmp1, dirDataTmp2, fileList, true);
    timer.toc("     alignment time: ");

    stopOpenMP();

    // ---------------------
    // Save offsets
    // ---------------------
    cResultAligner<float>().writeResultRigid(corrRes, dirDataOut + "/alignData.txt");
}

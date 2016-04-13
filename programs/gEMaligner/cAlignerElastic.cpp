/***********************************************************************
 *  File:       cAlignerElastic.cpp
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

void cAligner::alignElastic(void)
{
    // output directory
    dirDataTmp1 = dirDataOut + "/align";
    dirDataTmp2 = dirDataOut + "/tps";

    cFileSystem     objFileSystem(dirDataTmp2);

    if (objFileSystem.isExist()) {
        objFileSystem.dirDelete();
    }
    objFileSystem.dirCreate();

    // variables
    cTime                   timer;

    startOpenMP();

// tracking time: 2797 (s)
// tracking time without correlation: 19 (s)
// matrix time Armadillo: 755 (s)
// matrix time Eigen:     227 (s)
// matrix time without inversion: 14 (s)
// warping time GPU: 1913 (s)
// warping time 32 CPU: 31121

    // ---------------------
    // Point tracking
    // ---------------------
/*
    cSize3                  tplSizeHalf(10,10,0),
                            searchRange(10,10,0),
                            step(10,10,0),
                            border = tplSizeHalf + searchRange + 1;
*/
    cSize3                  tplSizeHalf(50,50,0),
                            searchRange(30,30,0),
                            step(50,50,0),
                            border = tplSizeHalf + searchRange + 1;

    cStackRegTPS2<float>    objRegTPS2;

    objRegTPS2.setParamsPointTracking(tplSizeHalf, searchRange,
                                      step, border);

    mpiMasterPrint(mpiID, "Tracking points...\n");
    timer.tic();
    objRegTPS2.trackPointsUsingCorr(dirDataTmp1, dirDataTmp2, fileList);
    timer.toc("     tracking time: ");

    // ---------------------
    // TPS2
    // ---------------------
    mpiMasterPrint(mpiID, "Computing TPS2...\n");
    timer.tic();
    cStackRegTPS2<double>().computeMatrixW(dirDataTmp2, fileList);
    timer.toc("     matrix time: ");

    // ---------------------
    // Warping
    // ---------------------
    mpiMasterPrint(mpiID, "Computing mapping...\n");
    timer.tic();
    cStackRegTPS2<double>().computeMappingAndInterp(dirDataTmp1,
                                dirDataTmp2, dirDataOut, fileList);
    timer.toc("     warping time: ");

    stopOpenMP();
}

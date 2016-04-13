/***********************************************************************
 *  File:       main.cpp
 *
 *  Purpose:    Entry main function for gEMfitter
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cFitter.hpp"

int main(int argc, char** argv)
{
    cFitter     fitter(argc, argv);
    cTime       timer;

    // pre-processing
    if (fitter.getMpiID() == 0) timer.tic();

    fitter.argParse();
    fitter.argVerify();
    fitter.argPrint();

    fitter.dataReadInput();

    if (fitter.getMpiID() == 0) timer.toc("Pre-processing time: ");
    mpiMasterPrint(fitter.getMpiID(), "\n");

    // correlation
    fitter.stopOpenMP();

    if (fitter.getComputeValue()) {
    //if (1) {

        mpiMasterPrint(fitter.getMpiID(), "Correlation calculation using " + num2str(fitter.getNumOrientation()) + " orientations...\n");
        if (fitter.getMpiID() == 0) timer.tic();

        fitter.computeCorr();
        fitter.dataWriteCorr();

        if (fitter.getMpiID() == 0) timer.toc("\n\nCorrelation time: ");
        mpiMasterPrint(fitter.getMpiID(), "\n");
    }
    else {
        fitter.dataReadCorr();
    }

    fitter.startOpenMP();

    // post-processing
    if (fitter.getMpiID() == 0 && fitter.getNumComponent() > 0) {

        mpiMasterPrint(fitter.getMpiID(), "Placing PDB using correlation maps...\n\n");

        if (fitter.getMpiID() == 0) timer.tic();

        fitter.computeResult();

        if (fitter.getMpiID() == 0) timer.toc("Placement time: ");
        mpiMasterPrint(fitter.getMpiID(), "\n");

        if (fitter.getChainName() != "none") {

            mpiMasterPrint(fitter.getMpiID(), "Computing RMSD...\n");

            fitter.computeRMSD();
        }
    }

    // final steps
    fitter.memFree();
    fitter.checkMemoryLeak();

    mpiMasterPrint(fitter.getMpiID(), "gEMfitter runs successfully!!!\n\n");

    return 0;
}

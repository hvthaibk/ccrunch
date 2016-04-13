/***********************************************************************
 *  File:       main.cpp
 *
 *  Purpose:    Entry main function for gEMaligner
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2014 Thai V. Hoang, IGBMC
 **********************************************************************/

#include "cAligner.hpp"

/* double
CPU time: 256.534
CPU time: 11.7129
GPU time: 0.648316
XmapDiff.getMax() = 1.05865e-09
XmapDiff.getMin() = -1.08139e-09
timeCPU2/timeCPU1 = 0.0456581
timeCPU1/timeGPU = 395.693
timeCPU2/timeGPU = 18.0666
*/
/* float
CPU time: 221.364
CPU time: 8.62353
GPU time: 0.357493
XmapDiff.getMax() = 0.717773
XmapDiff.getMin() = -0.662109
timeCPU2/timeCPU1 = 0.0389564
timeCPU1/timeGPU = 619.212
timeCPU2/timeGPU = 24.1222
*/

int main(int argc, char** argv)
{
    cAligner    aligner(argc, argv);
    cTime       timer;

    // pre-processing
    aligner.argParse();
    aligner.argVerify();
    aligner.argPrint();

    // remove curtaining
    if (aligner.getMode() == 0) {
        aligner.removeCurtain();
    }

    // translational alignment
    if (aligner.getMode() == 1) {
        aligner.alignRigid();
    }

    // elastic alignment
    if (aligner.getMode() == 2) {
        aligner.alignElastic();
    }

    // final steps
    aligner.memFree();
    aligner.checkMemoryLeak();

    mpiMasterPrint(aligner.getMpiID(), "gEMaligner runs successfully!!!\n\n");

    return 0;
}

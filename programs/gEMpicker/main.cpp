/***********************************************************************
 *  File:       main.cpp
 *
 *  Purpose:    Entry main function for gEMpicker
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "gEMpicker.hpp"

static sPickerParams    paramsPicker;
static cSystem          mySystem;
static cTime            timer;

int main(int argc, char** argv)
{
#ifdef __GEM_USE_MPI__
    // MPI Initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &paramsPicker.mpiNumProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &paramsPicker.mpiID);
#else
    paramsPicker.mpiID = 0;
    paramsPicker.mpiNumProcs = 1;
#endif

    mpiMasterPrint(paramsPicker.mpiID, "\n------------------------------------------------------------------------\n");
    mpiMasterPrint(paramsPicker.mpiID, "                Automatic Particle Picking by gEMpicker                   \n");
    mpiMasterPrint(paramsPicker.mpiID, "------------------------------------------------------------------------\n\n");

    // params processing
    if (paramsPicker.mpiID == 0) timer.tic();
    paramsInput(paramsPicker, argc, argv);
    mySystem.checkComputingResource(paramsPicker.nGPU, paramsPicker.nCPU);

    paramsVerify(paramsPicker);

    writeListRef(paramsPicker);
    writeListTpl(paramsPicker);
    writeListMsk(paramsPicker);

    if (paramsPicker.mpiID == 0) {
        paramsPrint(paramsPicker, argv[0]);
        timer.toc("Params processing time: ");
        std::cout << std::endl;
    }

    if (paramsPicker.mode == 0) {
#ifdef __GEM_USE_OPENMP__       // required
        omp_set_num_threads(1);
#endif

        // correlation
        mpiMasterPrint(paramsPicker.mpiID, "Computing correlation maps...\n\n");
        if (paramsPicker.mpiID == 0) timer.tic();
        threadFuncInit(paramsPicker);
        if (paramsPicker.mpiID == 0) timer.toc("Total computation time: ");
        mpiMasterPrint(paramsPicker.mpiID, "\n");
    }
    else {
#ifdef __GEM_USE_OPENMP__
        omp_set_num_threads(omp_get_num_procs());
#endif

        // picking
        mpiMasterPrint(paramsPicker.mpiID, "Picking particles using correlation maps...\n\n");
        if (paramsPicker.mpiID == 0) {
            timer.tic();
            pickParticle(paramsPicker);
            timer.toc("Picking time: ");
            mpiMasterPrint(paramsPicker.mpiID, "\n");
        }
    }

    /*if (paramsPicker.mpiID == 0) {
        msgTotalTime(paramsPicker, timer);
    }*/

#ifdef __GEM_USE_MPI__
    // MPI termination
    MPI_Finalize();
#endif

    mySystem.checkMemoryLeak();

    mpiMasterPrint(paramsPicker.mpiID, "gEMpicker runs successfully!!!\n\n");

    return 0;
}

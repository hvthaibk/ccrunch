/***********************************************************************
 *  File:       cFitterMem.cpp
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

cFitter::cFitter(int& argc, char**& argv)
{
    myArgc = argc;
    myArgv = new char*[myArgc];
    for(int i = 0; i < myArgc; i++) {
        myArgv[i] = new char[strlen(argv[i])+1];
        array_memcpy(myArgv[i], argv[i], strlen(argv[i])+1);
    }

#ifdef __GEM_USE_MPI__
    // MPI Initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiNumProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiID);
#else
    mpiID = 0;
    mpiNumProcs = 1;
#endif

    mpiMasterPrint(mpiID, "\n------------------------------------------------------------------------\n");
    mpiMasterPrint(mpiID, "            Rigid-Body Multiresolution Fitting by gEMfitter               \n");
    mpiMasterPrint(mpiID, "------------------------------------------------------------------------\n\n");
}

cFitter::~cFitter()
{
    for(int i = 0; i < myArgc; i++) {
        delete myArgv[i];
    }
    delete [] myArgv;

#ifdef __GEM_USE_MPI__
    // MPI termination
    MPI_Finalize();
#endif
}

void cFitter::memFree(void)
{
    mapRef.memFree();
    mapTpl.memFree();
    mapCor.memFree();
    mapInd.memFree();

    mapRefSurface.memFree();
    mapRefOuter.memFree();

    mapTplSurface.memFree();
    mapTplInner.memFree();
}

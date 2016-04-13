/***********************************************************************
 *  File:       cAlignerMem.cpp
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

cAligner::cAligner(int& argc, char**& argv)
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
    mpiMasterPrint(mpiID, "            Elastic 3D FIB/SEM Reconstruction by gEMaligner               \n");
    mpiMasterPrint(mpiID, "------------------------------------------------------------------------\n\n");
}

cAligner::~cAligner()
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

void cAligner::memFree(void)
{
/*    mapRef.memFree();
    mapTpl.memFree();
    mapCor.memFree();
    mapInd.memFree();

    mapRefSurface.memFree();
    mapRefOuter.memFree();

    mapTplSurface.memFree();
    mapTplInner.memFree();*/
}

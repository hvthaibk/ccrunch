/***********************************************************************
 *  File:       cAligner.hpp
 *
 *  Purpose:    Header file for cAligner
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2014 Thai V. Hoang, IGBMC
 **********************************************************************/

#ifndef __GEM_CFITTER_HPP__
#define __GEM_CFITTER_HPP__

#include "cThreadBoost.hpp"

#ifdef __GEM_USE_OPENMP__
#include <omp.h>
#endif

#ifdef __GEM_USE_MPI__
#include "mpi_xcorr.hpp"
#endif

#ifdef __GEM_USE_CUDA__
#include "cuCCorr.cuh"
#include "cuECorr.cuh"
#include "cuNCorr.cuh"
#include "cuWCorr.cuh"
#include "cuXCorr.cuh"
#include "cuData3.hpp"
#include "cuTime.hpp"
#include "cuRegTPS2.hpp"

#include "cuFilterFrequency.hpp"
#include "cuFilterSpatial.hpp"
#endif

#include "cCCorr.hpp"
#include "cECorr.hpp"
#include "cNCorr.hpp"
#include "cWCorr.hpp"
#include "cXCorr.hpp"
#include "cData3XCorr.hpp"
#include "cSystem.hpp"
#include "cTime.hpp"
#include "cString.hpp"
#include "cSetPoint3.hpp"
#include "cResultAligner.hpp"
#include "cFileSystem.hpp"

#include "cCCP4.hpp"
#include "cMRC.hpp"
#include "cTIFF.hpp"
#include "cOpenCV.hpp"

#include "cRegTPS2.hpp"
#include "cRegTrackPoint2.hpp"
#include "cStackRegTrans.hpp"
#include "cStackRegTPS2.hpp"

#include "cFilterFrequency.hpp"
#include "cFilterSpatial.hpp"

#include "cTransSimilarity.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

#include "tclap/CmdLine.h"

using namespace gem;


class cAligner : public cSystem
{
    int             mpiID, mpiNumProcs;
    int             nGPU,  nCPU;
    std::string     dirDataIn;
    std::string     dirDataOut;

    // internal variables
    int                           myArgc;
    char                        **myArgv;

    std::string                 dirDataTmp1, dirDataTmp2;
    std::vector<std::string>    fileList;
    unsigned int                mode;

public:
     cAligner(int& argc, char**& argv);
    ~cAligner();

    int             getMpiID (void) const { return mpiID; }
    unsigned int    getMode  (void) const { return mode;  }

    void    memFree  (void);

    void    argParse (void);
    void    argVerify(void);
    void    argPrint (void);

    void    removeCurtain(void);
    void    alignRigid   (void);
    void    alignElastic (void);
};

#endif

/***********************************************************************
 *  File:       gEMpicker.hpp
 *
 *  Purpose:    Header file for gEMpicker
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cThreadBoost.hpp"

#ifdef __GEM_USE_OPENMP__
#include <omp.h>
#endif

#ifdef __GEM_USE_MPI__
#include "mpi_xcorr.hpp"
#endif

#ifdef __GEM_USE_CUDA__
#include "cuNCorr.cuh"
#endif

#include "cNCorr.hpp"
#include "cData3XCorr.hpp"
#include "cFileSystem.hpp"
#include "cSystem.hpp"
#include "cTime.hpp"
#include "cString.hpp"
#include "cResultPicker.hpp"
#include "cTIFF.hpp"
#include "cOpenCV.hpp"
#include "cMRC.hpp"
#include "cTransSimilarity.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

#include "tclap/CmdLine.h"

using namespace gem;

struct sPickerParams
{
    std::string                 dirRefStr,    dirTplStr,
                                dirMskRefStr, dirMskTplStr,
                                dirResStr;
    unsigned int                mode;
    float                       angle2D;
    bool                        contrast;
    float                       threshold;
    float                       thresholdHigh;
    unsigned int                nPickMax;
    size_t                      boxSize;
    size_t                      boxDist;
    size_t                      boxBorder;
    int                         nCPU, nGPU;
    bool                        mpiDistTgt;
    int                         mpiID, mpiNumProcs;

    // internal variables
    bool                        bMskSingle;
    size_t                      numRot2D;
    std::vector<std::string>    fileListRef, fileListTpl, fileListMsk;
};

void    paramsInput (sPickerParams& params, int& argc, char**& argv);
void    paramsVerify(sPickerParams& params);
void    paramsPrint (sPickerParams& params, char* argv);
void    writeListRef(sPickerParams& params);
void    writeListTpl(sPickerParams& params);
void    writeListMsk(sPickerParams& params);
void    readListRef (sPickerParams& params);
void    readListTpl (sPickerParams& params);
void    readListMsk (sPickerParams& params);
void    msgTotalTime(sPickerParams& params, cTime& cTimer);

void    threadFuncInit(sPickerParams& params);

void    pickParticle  (sPickerParams& params);

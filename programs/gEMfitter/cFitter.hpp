/***********************************************************************
 *  File:       cFitter.hpp
 *
 *  Purpose:    Header file for cFitter
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
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
#include "cuTransThreshold.hpp"
#endif

#include "cCCorr.hpp"
#include "cECorr.hpp"
#include "cNCorr.hpp"
#include "cWCorr.hpp"
#include "cXCorr.hpp"
#include "cData3XCorr.hpp"
#include "cFileSystem.hpp"
#include "cSystem.hpp"
#include "cTime.hpp"
#include "cString.hpp"
#include "cPDB.hpp"
#include "cESBTL.hpp"
#include "cCCP4.hpp"
#include "cRandom.hpp"
#include "cSetPoint3.hpp"
#include "cResultFitter.hpp"
#include "cFilterSpatial.hpp"

#include "cTransDistance.hpp"
#include "cTransMorpho.hpp"
#include "cTransSimilarity.hpp"
#include "cTransThreshold.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

#include "pft_vertices.hpp"

#include "tclap/CmdLine.h"

using namespace gem;

extern const char *symTypeList[];

class cFitter : public cSystem
{
    int             mpiID, mpiNumProcs;
    int             nGPU,  nCPU;
    size_t          nNode, nRotZ;
    bool            bRotTexture;
    bool            bLaplacian;
    bool            bNyquist;
    double          resoRef;
    unsigned int    mCorr;
    unsigned int    mWCC;
    double          wDensity, wSurface, wPenalty;
    double          threshMsk;
    std::string     fileNameRef, fileNameTpl, fileNameMsk, dirNameRes;
    std::string     symType;
    bool            bRefine;
    double          tolRefine;
    unsigned int    iterRefine;
    unsigned int    nComponent;
    std::string     fileNameMdl;
    std::string     chainName;

    // internal variables
    int                       myArgc;
    char                    **myArgv;

    std::string             fileNameResCorr, fileNameResIndx, dirNameFit;

    cPDB                    pdbTpl, pdbMdl;
    cMapHeader              mapRefHeader, mapTplHeader, mapMskHeader;
    cData3XCorr<float>      mapRef, mapTpl, mapCor, mapInd;
    cVector3<size_t>        sizeRefOriginal;
    uint32_t                compute;

    // combined correlation
    cData3XCorr<float>      mapRefSurface, mapRefOuter;
    cData3XCorr<float>      mapTplSurface, mapTplInner;

public:
     cFitter(int& argc, char**& argv);
    ~cFitter();

    int             getMpiID         (void) const { return mpiID;       }
    uint32_t        getComputeValue  (void) const { return compute;     }
    unsigned int    getNumComponent  (void) const { return nComponent;  }
    size_t          getNumOrientation(void) const { return nNode*nRotZ; }
    std::string     getChainName     (void) const { return chainName;   }

    void    memFree(void);

    void    argParse (void);
    void    argVerify(void);
    void    argPrint (void);

    void    dataReadInput(void);
    void    dataReadCorr (void);
    void    dataWriteCorr(void);
    void    dataWriteTmp (void);

    void    computeCorr  (void);
    void    computeResult(void);
    void    computeRMSD  (void);

    double  refineNR3(cVector3<double>& transDouble, cVector3<double>& angleDouble);
    double  refineGSL(cVector3<double>& transDouble, cVector3<double>& angleDouble);
};

#endif

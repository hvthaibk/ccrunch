/***********************************************************************
 *  File:       cAlignerArg.cpp
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

void cAligner::argParse(void)
{
    const std::string   funcName("void cAligner::argParse(void)");

    try {
        TCLAP::CmdLine cmd("gEMaligner is an elastic 3D recontruction of FIB/SEM images. \n "
                           "Supported formats: TIFF \n "
                           "Copyright: Thai V. Hoang, IGBMC \n "
                           "Email: hvthaibk@gmail.com \n "
                           "Reference: \n Manuscript in preparation.", '=', "0.9");

        TCLAP::ValueArg<std::string>    arg_dirDataIn ("", "dirDataIn",  "Directory of input stack",                         true,    "", "string");
        TCLAP::ValueArg<std::string>    arg_dirDataOut("", "dirDataOut", "Directory of output stack",                        true,    "", "string");
        TCLAP::ValueArg<unsigned int>   arg_mode      ("", "mode",       "Calculation mode: 0-curtain, 1-rigid, 2-elastic",  true,     1, "unsigned int");
        TCLAP::ValueArg<unsigned int>   arg_nCPU      ("", "nCPU",       "Number of CPU cores / process: (1)",              false,     1, "unsigned int");
        TCLAP::ValueArg<unsigned int>   arg_nGPU      ("", "nGPU",       "Number of GPUs / process: (0)",                   false,     0, "unsigned int");

        cmd.add(arg_nGPU);
        cmd.add(arg_nCPU);
        cmd.add(arg_mode);
        cmd.add(arg_dirDataOut);
        cmd.add(arg_dirDataIn);

        cmd.parse(myArgc, myArgv);

        dirDataIn  = arg_dirDataIn.getValue();
        dirDataOut = arg_dirDataOut.getValue();
        mode       = arg_mode.getValue();
        nCPU       = arg_nCPU.getValue();
        nGPU       = arg_nGPU.getValue();

        checkComputingResource(nGPU, nCPU);
    }
    catch (TCLAP::ArgException &e) {
        if (mpiID == 0) {
            std::cerr << "Parsing error in " + funcName + ": "
                      << e.error()
                      << " for arg "
                      << e.argId() << std::endl;
        }

        exit(EXIT_FAILURE);
    }
}

void cAligner::argVerify(void)
{
    const std::string   funcName("void cAligner::argVerify(void)");

    // input stack

    require(cFileSystem(dirDataIn).isExist(),
            funcName + ": " + dirDataIn + " does not exist");

    cFileSystem(dirDataIn).getFileList(".tif", fileList);
    require(fileList.size() > 0, funcName + ": input image stack is empty");

    // output directory
    cFileSystem(dirDataOut).dirCreate();
}

void cAligner::argPrint(void)
{
    const std::string   funcName("void cAligner::argPrint(void)");

    if (mpiID != 0) return;

    std::string     strGap = "\n               ";

    std::cout << "Running " << myArgv[0] << ": use ";
#ifdef __GEM_USE_MPI__
    std::cout << "MPI, ";
#endif
#ifdef __GEM_USE_OPENMP__
    std::cout << "OpenMP, ";
#endif
#ifdef __GEM_USE_CUDA__
    std::cout << "CUDA, ";
#endif

    std::cout << strGap << "number of processes (MPI)     = " << mpiNumProcs
              << strGap << "number of GPUs      / process = " << nGPU
              << strGap << "number of CPU cores / process = " << nCPU;

    std::cout << strGap << "input directory               = " << dirDataIn;
    std::cout << strGap << "output directory              = " << dirDataOut;
    std::cout << strGap << "calculation mode              = " << mode;
    switch (mode) {
        case 0:     std::cout << " (curtain removal)";       break;
        case 1:     std::cout << " (rigid alignment)";       break;
        case 2:     std::cout << " (elastic registration)";  break;
        default:    ERROR(funcName, "unsupported calculation mode");
    }

    std::cout << std::endl << std::endl;
}

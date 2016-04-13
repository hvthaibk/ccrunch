/***********************************************************************
 *  File:       cFitterArg.cpp
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

const char *symTypeList[] = {"C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
                             "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10"};

template<typename T, size_t N>
T* end(T (&ra)[N]) { return ra + N; }

void cFitter::argParse(void)
{
    const std::string   funcName("void cFitter::argParse(void)");

    try {
        TCLAP::CmdLine cmd("gEMfitter is a template-based, highly parallel GPU-accelerated program for "
                           "multi-resolution fitting of macromolecular structures developed at INRIA Nancy. "
                           "It can run in multi-thread and multi-process mode to use all available CPU cores and GPUs "
                           "in a single workstation or a multi-node cluster "
                           "in order to speedup the calculation of the scoring functions. \n "
                           "Supported formats: CCP4 (target, mask, result) and PDB (search, model) \n "
                           "Copyright: Thai V. Hoang, INRIA \n "
                           "Email: hvthaibk@gmail.com \n "
                           "Reference: \n T. V. Hoang, X. Cavin, and D. W. Ritchie, "
                           "\"gEMfitter: A highly parallel FFT-based 3D density fitting tool with GPU texture memory acceleration,\" "
                           "Journal of Structural Biology, 2013.", '=', "1.1");

        TCLAP::ValueArg<size_t>         arg_nVertice   ("", "nVertice",    "Number of icosahedral vertices",                             true,       0, "unsigned int");
        TCLAP::ValueArg<bool>           arg_bNyquist   ("", "bNyquist",    "Map resampling using the Nyquist rate: (0)-no, 1-yes",       false,      0, "bool");
        TCLAP::ValueArg<bool>           arg_bLaplacian ("", "bLaplacian",  "Laplacian operation: (0)-no, 1-yes",                         false,      0, "bool");
        TCLAP::ValueArg<size_t>         arg_nRotZ      ("", "nRotZ",       "Number of z-axis rotations",                                 true,       0, "unsigned int");
        TCLAP::ValueArg<bool>           arg_bRotTexture("", "bRotTexture", "GPU 3D rotation using texture fetching: 0-no, (1)-yes",      false,      1, "bool");
        TCLAP::ValueArg<unsigned int>   arg_mCorr      ("", "mCorr",       "Correlation mode: 0-XCC, 1-ECC, 2-NCC, 3-WCC, 4-SCC, 5-CCC", true,       0, "unsigned int");
        TCLAP::ValueArg<unsigned int>   arg_mWCC       ("", "mWCC",        "WCC definition: 0-constant NCC, (1)-distance",               false,      1, "unsigned int");
        TCLAP::ValueArg<double>         arg_wDensity   ("", "wDensity",    "Weight for density correlation in CCC mode (0.5)",           false,    0.5, "double");
        TCLAP::ValueArg<double>         arg_wSurface   ("", "wSurface",    "Weight for surface correlation in CCC mode (0.5)",           false,    0.5, "double");
        TCLAP::ValueArg<double>         arg_wPenalty   ("", "wPenalty",    "Weight for penalty correlation in CCC mode (0.5)",           false,    0.5, "double");
        TCLAP::ValueArg<std::string>    arg_fileRef    ("", "fileTgt",     "Target map",                                                 true,      "", "string");
        TCLAP::ValueArg<std::string>    arg_fileTpl    ("", "fileSch",     "Search pdb",                                                 true,      "", "string");
        TCLAP::ValueArg<std::string>    arg_fileMsk    ("", "fileMsk",     "Mask map (\"\")",                                            false, "none", "string");
        TCLAP::ValueArg<std::string>    arg_dirRes     ("", "dirRes",      "Result directory",                                           true,      "", "string");
        TCLAP::ValueArg<double>         arg_resoRef    ("", "resoTgt",     "Target map's resolution",                                    true,       0, "double");
        TCLAP::ValueArg<double>         arg_threshMsk  ("", "threshMsk",   "Thresholding value for creating mask (-10)",                 false,    -10, "double");
        TCLAP::ValueArg<std::string>    arg_symType    ("", "symType",     "Symmetry type (\"\")",                                       false,     "", "string");
        TCLAP::ValueArg<bool>           arg_bRefine    ("", "bRefine",     "Off-latice refinement: (0)-no, 1-yes",                       false,      0, "bool");
        TCLAP::ValueArg<double>         arg_tolRefine  ("", "tolRefine",   "Off-latice refinement: tolerance (1e-6)",                    false,   1e-6, "double");
        TCLAP::ValueArg<unsigned int>   arg_iterRefine ("", "iterRefine",  "Off-latice refinement: #iteration (25)",                     false,     25, "unsigned int");
        TCLAP::ValueArg<unsigned int>   arg_nComponent ("", "nComponent",  "Number of components: (0)",                                  false,      0, "unsigned int");
        TCLAP::ValueArg<std::string>    arg_fileMdl    ("", "fileMdl",     "Model PDB (\"\")",                                           false, "none", "string");
        TCLAP::ValueArg<std::string>    arg_chainName  ("", "chainName",   "Comparing chains (\"\")",                                    false, "none", "string");
        TCLAP::ValueArg<unsigned int>   arg_nCPU       ("", "nCPU",        "Number of CPU cores / process: (1)",                         false,      1, "unsigned int");
        TCLAP::ValueArg<unsigned int>   arg_nGPU       ("", "nGPU",        "Number of GPUs / process: (0)",                              false,      0, "unsigned int");

        cmd.add(arg_nGPU);
        cmd.add(arg_nCPU);
        cmd.add(arg_chainName);
        cmd.add(arg_fileMdl);
        cmd.add(arg_nComponent);
        cmd.add(arg_iterRefine);
        cmd.add(arg_tolRefine);
        cmd.add(arg_bRefine);
        cmd.add(arg_symType);
        cmd.add(arg_threshMsk);
        cmd.add(arg_resoRef);
        cmd.add(arg_dirRes);
        cmd.add(arg_fileMsk);
        cmd.add(arg_fileTpl);
        cmd.add(arg_fileRef);
        cmd.add(arg_wPenalty);
        cmd.add(arg_wSurface);
        cmd.add(arg_wDensity);
        cmd.add(arg_mWCC);
        cmd.add(arg_mCorr);
        cmd.add(arg_bRotTexture);
        cmd.add(arg_bLaplacian);
        cmd.add(arg_bNyquist);
        cmd.add(arg_nRotZ);
        cmd.add(arg_nVertice);

        cmd.parse(myArgc, myArgv);

        nNode       = arg_nVertice.getValue();
        nRotZ       = arg_nRotZ.getValue();
        bNyquist    = arg_bNyquist.getValue();
        bLaplacian  = arg_bLaplacian.getValue();
        bRotTexture = arg_bRotTexture.getValue();
        mCorr       = arg_mCorr.getValue();
        mWCC        = arg_mWCC.getValue();
        wDensity    = arg_wDensity.getValue();
        wSurface    = arg_wSurface.getValue();
        wPenalty    = arg_wPenalty.getValue();
        fileNameRef = arg_fileRef.getValue();
        fileNameTpl = arg_fileTpl.getValue();
        fileNameMsk = arg_fileMsk.getValue();
        dirNameRes  = arg_dirRes.getValue();
        resoRef     = arg_resoRef.getValue();
        threshMsk   = arg_threshMsk.getValue();
        symType     = arg_symType.getValue();
        bRefine     = arg_bRefine.getValue();
        tolRefine   = arg_tolRefine.getValue();
        iterRefine  = arg_iterRefine.getValue();
        nComponent  = arg_nComponent.getValue();
        fileNameMdl = arg_fileMdl.getValue();
        chainName   = arg_chainName.getValue();
        nCPU        = arg_nCPU.getValue();
        nGPU        = arg_nGPU.getValue();

        if (dirNameRes[dirNameRes.length()-1] != '/') {
            dirNameRes += "/";
        }

        if (mCorr == CC_MODE_SCORR) {
            require(bLaplacian == 0, funcName + ": surface fitting (SCC) cannot be use with Laplacian");
        }

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

void cFitter::argVerify(void)
{
    const std::string   funcName("void cFitter::argVerify(void)");

    std::string     message;
    std::string     hostname = getComputerName();
    std::string     fileNameRes;

    // reference
    message = "(" + hostname + ") ";
    message += "target map " + fileNameRef + " does not exist";
    require(cFileSystem(fileNameRef).isExist(), message.c_str());
    mpiMasterPrint(mpiID, "Target: " + fileNameRef + "\n");

    // template
    message = "(" + hostname + ") ";
    message += "search map " + fileNameTpl + " does not exist";
    require(cFileSystem(fileNameTpl).isExist(), message.c_str());
    mpiMasterPrint(mpiID, "Search: " + fileNameTpl + "\n");

    // mask
    if (fileNameMsk != "none") {
        message = "(" + hostname + ") ";
        message += "mask map " + fileNameMsk + " does not exist";
        require(cFileSystem(fileNameMsk).isExist(), message.c_str());
        mpiMasterPrint(mpiID, "Mask:   " + fileNameMsk + "\n");
    }

    // PDB model
    if (fileNameMdl != "none") {
        message = "(" + hostname + ") ";
        message += "pdb model " + fileNameMdl + " does not exist";
        require(cFileSystem(fileNameMdl).isExist(), message.c_str());
        mpiMasterPrint(mpiID, "Model:  " + fileNameMdl + "\n");
    }

    // result
    if (dirNameRes[dirNameRes.length()-1] != '/') {
        dirNameRes += "/";
    }

    fileNameRes = cFileSystem(fileNameRef).getFileRoot() +
           "__" + cFileSystem(fileNameTpl).getFileRoot() +
           "__" + num2str(nNode) + "__" + num2str(nRotZ);

    switch (mCorr) {
        case CC_MODE_XCORR: fileNameRes += "__XCC";                         break;
        case CC_MODE_ECORR: fileNameRes += "__ECC";                         break;
        case CC_MODE_NCORR: fileNameRes += "__NCC";                         break;
        case CC_MODE_WCORR: fileNameRes += "__WCC"   + num2str(mWCC);       break;
        case CC_MODE_SCORR: fileNameRes += "__SCC";                         break;
        case CC_MODE_CCORR: fileNameRes += "__CCC"   + num2str(mWCC)
                                              + "__" + num2str(wDensity,1)
                                              + "__" + num2str(wSurface,1)
                                              + "__" + num2str(wPenalty,1); break;
        default:            ERROR(funcName, "unsupported correlation mode");
    }

    fileNameRes += "__" + num2str(bLaplacian);
    if (fileNameMsk == "none") {
        fileNameRes += "__" + num2str(threshMsk,4);
    }
    else {
        fileNameRes += "__file";
    }
    mpiMasterPrint(mpiID, "Result: " + dirNameRes + fileNameRes + "\n");

    dirNameFit      = dirNameRes + fileNameRes + "/";
    fileNameResCorr = dirNameFit + "resCorr.map";
    fileNameResIndx = dirNameFit + "resIndx.map";

    if (!cFileSystem(dirNameRes).isExist()) {
        cFileSystem(dirNameRes).dirCreate(hostname);
    }
    if (!cFileSystem(dirNameFit).isExist()) {
        cFileSystem(dirNameFit).dirCreate(hostname);
    }

    // symmetry
    std::vector<std::string>    symVec(symTypeList, end(symTypeList));
    for (size_t i = 0; i < symVec.size(); i++) {
        if (symType == symVec[i]) {
        }
    }

    // resolution
    message = "(" + hostname + ") ";
    message += "target resolution must be in the range [2 50]";
    require(resoRef >= 2 && resoRef <= 50, message.c_str());

    // other options
    pft_vertices_check(nNode);

    // compute correlation?
    mpiMasterPrint(mpiID, "\n");
    if (mpiID == 0 && cFileSystem(fileNameResCorr).isExist()) {
        message = "(" + hostname + ") ";
        message += fileNameResCorr + " exists, no need to compute correlation maps again!!!\n\n";
        std:: cout << message;
        compute = 0;
    }
    else {
        compute = 1;
    }

#ifdef __GEM_USE_MPI__
    MPI_Bcast(&compute, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#endif  // __GEM_USE_MPI__
}

void cFitter::argPrint(void)
{
    const std::string   funcName("void cFitter::argPrint(void)");

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
              << strGap << "number of CPU cores / process = " << nCPU
              << strGap << "number of icosahedral nodes   = " << nNode << " (" << num2str(pft_vertices_getangle(nNode),1) << " degree)"
              << strGap << "number of Z-rotations         = " << nRotZ << " (" << num2str(360/(float)nRotZ,1) << " degree)";

    switch (mCorr) {
        case CC_MODE_XCORR: std::cout << strGap << "correlation mode              = XCC";         break;
        case CC_MODE_ECORR: std::cout << strGap << "correlation mode              = ECC";         break;
        case CC_MODE_NCORR: std::cout << strGap << "correlation mode              = NCC";         break;
        case CC_MODE_WCORR: std::cout << strGap << "correlation mode              = WCC"
                                                << " (mWCC=" << mWCC << ")";                      break;
        case CC_MODE_SCORR: std::cout << strGap << "correlation mode              = SCC";         break;
        case CC_MODE_CCORR: std::cout << strGap << "correlation mode              = CCC"
                                                << " (wDensity=" << wDensity
                                                << ", wSurface=" << wSurface
                                                << ", wPenalty=" << wPenalty << ")";              break;
        default:            ERROR(funcName, "unsupported correlation mode");
    }

    if (bNyquist)    std::cout << strGap << "Nyquist-rate resampling       = Yes";
    else             std::cout << strGap << "Nyquist-rate resampling       = No ";
    if (bLaplacian ) std::cout << strGap << "Laplacian pre-filtering       = Yes";
    else             std::cout << strGap << "Laplacian pre-filtering       = No ";
    if (bRotTexture) std::cout << strGap << "texture rotation in GPU       = Yes";
    else             std::cout << strGap << "texture rotation in GPU       = No ";

    std::cout << strGap << "target resolution             = " << resoRef;
    std::cout << strGap << "mask thresholding value       = " << threshMsk;
    if (fileNameMsk != "none") {
        std::cout << " (no need since a mask map is provided)";
    }
    if (symType != "") std::cout << strGap << "symmetry type                 = " << symType;
    else               std::cout << strGap << "symmetry type                 = No";

    if (bRefine) {
        std::cout << strGap << "off-latice refinement         = Yes";
        std::cout << strGap << "correlation peak tolerance    = " << tolRefine;
        std::cout << strGap << "max number of iterations      = " << iterRefine;
    }
    else {
        std::cout << strGap << "off-latice refinement         = No ";
    }

    std::cout << strGap << "number of components          = " << nComponent;
    std::cout << strGap << "chain name                    = " << chainName;

    std::cout << std::endl << std::endl;
}

/***********************************************************************
 *  File:       processArgsList.cpp
 *
 *  Purpose:    Implementation of params-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "gEMpicker.hpp"

void paramsInput(sPickerParams& params, int& argc, char**& argv)
{
    try {
        TCLAP::CmdLine cmd("gEMpicker is a template-based, highly parallel GPU-accelerated particle picking program "
                           "developed at INRIA Nancy. It runs in multi-thread and multi-process mode to use "
                           "all available CPU cores and GPUs in a single workstation or a multi-node cluster "
                           "for the correlation computation. \n "
                           "Supported formats: MRC for micrographs and references, TIFF for masks \n "
                           "Copyright: Thai V. Hoang, INRIA \n "
                           "Email: hvthaibk@gmail.com \n "
                           "Reference: \n T. V. Hoang, X. Cavin, P. Schultz, and D. W. Ritchie, "
                           "\"gEMpicker: A highly parallel GPU-accelerated tool for electron micrograph particle picking,\" "
                           "BMC Structural Biology, 2013.", '=', "1.1");

        TCLAP::ValueArg<std::string>    arg_dirRef    ("", "dirTgt",     "Directory of target images / micrographs (*.mrc)",                                     true,  "", "string");
        TCLAP::ValueArg<std::string>    arg_dirTpl    ("", "dirSch",     "Directory of search images / references (*.mrc)",                                      true,  "", "string");
        TCLAP::ValueArg<std::string>    arg_dirMskRef ("", "dirMskTgt",  "Directory of masks for target images / micrographs (*.tif)",                           false, "", "string");
        TCLAP::ValueArg<std::string>    arg_dirMskTpl ("", "dirMskSch",  "Directory of masks for search images / references (*.tif)",                            true,  "", "string");
        TCLAP::ValueArg<std::string>    arg_dirRes    ("", "dirRes",     "Directory for results",                                                                true,  "", "string");
        TCLAP::ValueArg<unsigned int>   arg_mode      ("", "mode",       "Running mode: 0-compute correlation, 1-perform picking (box only), 2-perform picking", true,   0, "unsigned int");
        TCLAP::ValueArg<float>          arg_angle2D   ("", "angle2D",    "In-plane rotating angle in degree, (0 = no rotation)",                                 true,   0, "float");
        TCLAP::ValueArg<bool>           arg_contrast  ("", "contrast",   "Micrograph contrast: 0-negative peaks, 1-positive peaks",                              true,   1, "bool");
        TCLAP::ValueArg<float>          arg_thresh    ("", "thresh",     "Threshold value for picking - low limit: (0,1]",                                       true,   0, "float");
        TCLAP::ValueArg<float>          arg_threshHigh("", "threshHigh", "Threshold value for picking - high limit: (0,1] (default = 1)",                        false,  1, "float");
        TCLAP::ValueArg<unsigned int>   arg_nPickMax  ("", "nPickMax",   "Max number of particles picked from each micrograph (default = 0, no limit)",          false,  0, "unsigned int");
        TCLAP::ValueArg<unsigned int>   arg_boxSize   ("", "boxSize",    "Size of picked images (if not provided or = 0, use reference size)",                   false,  0, "unsigned int");
        TCLAP::ValueArg<unsigned int>   arg_boxDist   ("", "boxDist",    "Min distance between peaks (if not provided or = 0, use reference size)",              false,  0, "unsigned int");
        TCLAP::ValueArg<unsigned int>   arg_boxBorder ("", "boxBorder",  "Min distance from box to micrograph boundary (default = 0)",                           false,  0, "unsigned int");
        TCLAP::ValueArg<unsigned int>   arg_nCPU      ("", "nCPU",       "Number of CPU cores / process: (1)",                                                   false,  1, "unsigned int");
        TCLAP::ValueArg<unsigned int>   arg_nGPU      ("", "nGPU",       "Number of GPUs / process: (0)",                                                        false,  0, "unsigned int");
        TCLAP::ValueArg<bool>           arg_mpiDistTgt("", "mpiDistTgt", "Distribution of micrographs to processes: (0)-no, 1-yes",                              false,  0, "bool");

        cmd.add(arg_mpiDistTgt);
        cmd.add(arg_nGPU);
        cmd.add(arg_nCPU);
        cmd.add(arg_boxBorder);
        cmd.add(arg_boxDist);
        cmd.add(arg_boxSize);
        cmd.add(arg_nPickMax);
        cmd.add(arg_threshHigh);
        cmd.add(arg_thresh);
        cmd.add(arg_contrast);
        cmd.add(arg_angle2D);
        cmd.add(arg_mode);
        cmd.add(arg_dirRes);
        cmd.add(arg_dirMskTpl);
        cmd.add(arg_dirMskRef);
        cmd.add(arg_dirTpl);
        cmd.add(arg_dirRef);

        cmd.parse(argc, argv);

        params.dirRefStr     = arg_dirRef.getValue();
        params.dirTplStr     = arg_dirTpl.getValue();
        params.dirMskRefStr  = arg_dirMskRef.getValue();
        params.dirMskTplStr  = arg_dirMskTpl.getValue();
        params.dirResStr     = arg_dirRes.getValue();
        params.mode          = arg_mode.getValue();
        params.angle2D       = arg_angle2D.getValue();
        params.contrast      = arg_contrast.getValue();
        params.threshold     = arg_thresh.getValue();
        params.thresholdHigh = arg_threshHigh.getValue();
        params.nPickMax      = arg_nPickMax.getValue();
        params.boxSize       = arg_boxSize.getValue();
        params.boxDist       = arg_boxDist.getValue();
        params.boxBorder     = arg_boxBorder.getValue();
        params.nCPU          = arg_nCPU.getValue();
        params.nGPU          = arg_nGPU.getValue();
        params.mpiDistTgt    = arg_mpiDistTgt.getValue();

        if (params.dirRefStr[params.dirRefStr.length()-1] != '/')   params.dirRefStr += "/";
        if (params.dirTplStr[params.dirTplStr.length()-1] != '/')   params.dirTplStr += "/";
        if (params.dirResStr[params.dirResStr.length()-1] != '/')   params.dirResStr += "/";

        if (params.dirMskRefStr[params.dirMskRefStr.length()-1] != '/')   params.dirMskRefStr += "/";
        if (params.dirMskTplStr[params.dirMskTplStr.length()-1] != '/')   params.dirMskTplStr += "/";
    }
    catch (TCLAP::ArgException &e) {
        if (params.mpiID == 0) {
            std::cerr << "Parsing error in void paramsInput(sPickerParams& params, int& argc, char**& argv): "
                      << e.error()
                      << " for arg "
                      << e.argId() << std::endl;
        }

        exit(EXIT_FAILURE);
    }
}

void paramsVerify(sPickerParams& params)
{
    cSystem         sysInfo;
    std::string     message;
    std::string     hostname = sysInfo.getComputerName();
    std::string     dirName;

    // reference
    message = "(" + hostname + ") ";
    message += "micrograph directory " + params.dirRefStr + " does not exist";
    require(cFileSystem(params.dirRefStr).isExist(), message.c_str());

    cFileSystem(params.dirRefStr).getFileList(".mrc", params.fileListRef, true);

    message  = "(" + hostname + ") ";
    message += "micrograph directory " + params.dirRefStr + " contains no MRC file";
    require(params.fileListRef.size() > 0, message.c_str());

    mpiMasterPrint(params.mpiID, "Micrograph: " + params.dirRefStr + "\n");

    std::sort(params.fileListRef.begin(), params.fileListRef.end());

    // template
    message  = "(" + hostname + ") ";
    message += "reference directory " + params.dirTplStr + " does not exist";
    require(cFileSystem(params.dirTplStr).isExist(), message.c_str());

    cFileSystem(params.dirTplStr).getFileList(".mrc", params.fileListTpl, true);

    message  = "(" + hostname + ") ";
    message += "reference directory " + params.dirTplStr + " contains no MRC file";
    require(params.fileListTpl.size() > 0, message.c_str());

    mpiMasterPrint(params.mpiID, "Reference:  " + params.dirTplStr + "\n");

    std::sort(params.fileListTpl.begin(), params.fileListTpl.end());

    // template mask
    message  = "(" + hostname + ") ";
    message += "mask directory " + params.dirMskTplStr + " does not exist";
    require(cFileSystem(params.dirMskTplStr).isExist(), message.c_str());

    cFileSystem(params.dirMskTplStr).getFileList(".tif", params.fileListMsk, true);

    message  = "(" + hostname + ") ";
    message += "mask directory " + params.dirMskTplStr + " contains no TIFF file";
    require(params.fileListMsk.size() > 0, message.c_str());

    mpiMasterPrint(params.mpiID, "Mask:       " + params.dirMskTplStr + "\n");

    std::sort(params.fileListMsk.begin(), params.fileListMsk.end());

    // result
    dirName = params.dirResStr;
    if (!cFileSystem(dirName).isExist()) {
         cFileSystem(dirName).dirCreate(hostname);
    }

    dirName = params.dirResStr + "xcorr";
    if (!cFileSystem(dirName).isExist()) {
         cFileSystem(dirName).dirCreate(hostname);
    }

    dirName = params.dirResStr + "pik_coord";
    if (!cFileSystem(dirName).isExist()) {
         cFileSystem(dirName).dirCreate(hostname);
    }

    dirName = params.dirResStr + "pik_box";
    if (!cFileSystem(dirName).isExist()) {
         cFileSystem(dirName).dirCreate(hostname);
    }

    dirName = params.dirResStr + "pik_ext";
    if (!cFileSystem(dirName).isExist()) {
         cFileSystem(dirName).dirCreate(hostname);
    }

    mpiMasterPrint(params.mpiID, "Result:     " + params.dirResStr + "\n\n");

    // templates vs. masks
    std::string        fileNameTpl, fileNameMsk;
    if (params.fileListMsk.size() > 1) {
        message  = "(" + hostname + ") ";
        message += "#references != #mask images ";
        message += "(#Ref = "  + num2str(params.fileListTpl.size());
        message += ", #Msk = " + num2str(params.fileListMsk.size());
        message += ")";
        require(params.fileListTpl.size() == params.fileListMsk.size(), message.c_str());

        for (size_t i = 0; i < params.fileListTpl.size(); i++) {
            fileNameTpl = cFileSystem(params.fileListTpl[i]).getFileRoot();
            fileNameMsk = cFileSystem(params.fileListMsk[i]).getFileRoot();

            require(fileNameTpl == fileNameMsk, "reference fileName != mask fileName (" +
                                                fileNameTpl + " != " + fileNameMsk + ")");
        }
    }

    // other options
    if (params.mode) {    // picking mode
        require((params.threshold     > 0) && (std::abs(params.threshold)     <= 1) &&
                (params.thresholdHigh > 0) && (std::abs(params.thresholdHigh) <= 1),
                "thresholding value must be in (0 1]");
        require(params.threshold < params.thresholdHigh, "high threshold must be greater than low threshold");
    }
    require(params.angle2D >=  0, "in-plane rotating angle must be in [0 360)");
    require(params.angle2D < 360, "in-plane rotating angle must be in [0 360)");
    require((params.contrast == 0) || (params.contrast == 1),
            "contrast parameter must only take value 0 or 1");

    // internal variables
    if (params.fileListMsk.size() == 1) params.bMskSingle = true;
    else                                params.bMskSingle = false;

    params.numRot2D = (params.angle2D == 0) ? 1 : (size_t) std::floor(360/params.angle2D);
}

void paramsPrint(sPickerParams& params, char* argv)
{
    const std::string   funcName("void paramsPrint(sPickerParams& params, char* argv)");

    std::string     strGap = "\n               ";

    std::cout << "Running " << argv << ": use ";
#ifdef __GEM_USE_MPI__
    std::cout << "MPI, ";
#endif
#ifdef __GEM_USE_OPENMP__
    std::cout << "OpenMP, ";
#endif
#ifdef __GEM_USE_CUDA__
    std::cout << "CUDA, ";
#endif

    std::cout << strGap << "number of micrographs              = " << params.fileListRef.size()
              << strGap << "number of references               = " << params.fileListTpl.size()
              << strGap << "number of mask images              = " << params.fileListMsk.size();
    if (params.fileListMsk.size() == 1) {
        std::cout << " (identical mask)";
    }
    else {
        std::cout << " (individual masks)";
    }

    std::cout << strGap << "running mode                       = " << params.mode;
    if (params.mode == 0) {
        std::cout << " (correlation)";
    }
    else if (params.mode == 1) {
        std::cout << " (picking, box only)";
    }
    else if (params.mode == 2) {
        std::cout << " (picking)";
    }
    else {
        ERROR(funcName, "incorrect running mode");
    }

    std::cout << strGap << "in-plane rotating angle (degree)   = " << params.angle2D;
    if (params.angle2D == 0) {
        std::cout << " (no rotation)";
    }

    std::cout << strGap << "micrograph contrast                = " << params.contrast;
    if (params.contrast == 0) {
        std::cout << " (negative peaks)";
    }
    else {
        std::cout << " (positive peaks)";
    }

    std::cout << strGap << "peak threshold value (low)         = " << params.threshold;
    std::cout << strGap << "peak threshold value (high)        = " << params.thresholdHigh;

    std::cout << strGap << "max number of picked particles     = " << params.nPickMax;
    if (params.nPickMax == 0) {
        std::cout << " (no limit)";
        params.nPickMax = UINT_MAX;
    }

    std::cout << strGap << "size of picked particles           = " << params.boxSize;
    if (params.boxSize == 0) {
        std::cout << " (use reference size)";
    }

    std::cout << strGap << "min distance between peaks         = " << params.boxDist;
    if (params.boxDist == 0) {
        std::cout << " (use reference size)";
    }

    std::cout << strGap << "min distance box-micrograph        = " << params.boxBorder;

    std::cout << strGap << "number of processes (MPI)          = " << params.mpiNumProcs
              << strGap << "number of GPUs      / process      = " << params.nGPU
              << strGap << "number of CPU cores / process      = " << params.nCPU;

    if (params.mpiNumProcs > 1 && params.mpiDistTgt) {
        std::cout << strGap << "distribution of micrographs (MPI)  = Yes";
    }
    else {
        std::cout << strGap << "distribution of micrographs (MPI)  = No";
    }

    std::cout << strGap << std::endl;
}

void writeListRef(sPickerParams& params)
{
    std::ofstream       fileHandle;
    std::string         fileName = params.dirResStr + "ListTgt.txt";

    fileHandle.open(fileName.c_str());
    assure(fileHandle, "ListTgt.txt");

    fileHandle << "Total number of micrographs: " << params.fileListRef.size() << std::endl;

    for (size_t i = 0; i < params.fileListRef.size(); i++) {
        fileHandle << std::setw(5) << i+1 << "     "
                   << params.fileListRef[i]
                   << std::endl;
    }

    fileHandle.close();
}

void writeListTpl(sPickerParams& params)
{
    std::ofstream       fileHandle;
    std::string         fileName = params.dirResStr + "ListSch.txt";

    fileHandle.open(fileName.c_str());
    assure(fileHandle, "ListSch.txt");

    fileHandle << "Total number of references: " << params.fileListTpl.size() << std::endl;

    for (size_t i = 0; i < params.fileListTpl.size(); i++) {
        fileHandle << std::setw(5) << i+1 << "     "
                   << params.fileListTpl[i]
                   << std::endl;
    }

    fileHandle.close();
}

void writeListMsk(sPickerParams& params)
{
    std::ofstream       fileHandle;
    std::string         fileName = params.dirResStr + "ListMsk.txt";

    fileHandle.open(fileName.c_str());
    assure(fileHandle, "ListMsk.txt");

    fileHandle << "Total number of masks: " << params.fileListMsk.size() << std::endl;

    for (size_t i = 0; i < params.fileListMsk.size(); i++) {
        fileHandle << std::setw(5) << i+1 << "     "
                   << params.fileListMsk[i]
                   << std::endl;
    }

    fileHandle.close();
}

void readListRef(sPickerParams& params)
{
    size_t              sizeNum;
    std::string         line;
    std::ifstream       fileHandle;
    std::string         fileName;

    fileName = params.dirResStr + "ListTgt.txt";
    fileHandle.open(fileName.c_str());
    assure(fileHandle, "ListTgt.txt");

    getline(fileHandle, line);
    sizeNum = std::atoi(line.substr(23,line.length()).c_str());
    params.fileListRef.resize(sizeNum);

    for (size_t i = 0; i < sizeNum; i++) {
        getline(fileHandle, line);
        params.fileListRef[i] = line.substr(10,line.length());
    }

    fileHandle.close();
}

void readListTpl(sPickerParams& params)
{
    size_t              sizeNum;
    std::string         line;
    std::ifstream       fileHandle;
    std::string         fileName;

    fileName = params.dirResStr + "ListSch.txt";
    fileHandle.open(fileName.c_str());
    assure(fileHandle, "ListSch.txt");

    getline(fileHandle, line);
    sizeNum = std::atoi(line.substr(23,line.length()).c_str());
    params.fileListTpl.resize(sizeNum);

    for (size_t i = 0; i < sizeNum; i++) {
        getline(fileHandle, line);
        params.fileListTpl[i] = line.substr(10,line.length());
    }

    fileHandle.close();
}

void readListMsk(sPickerParams& params)
{
    size_t              sizeNum;
    std::string         line;
    std::ifstream       fileHandle;
    std::string         fileName;

    fileName = params.dirResStr + "ListMsk.txt";
    fileHandle.open(fileName.c_str());
    assure(fileHandle, "ListMsk.txt");

    getline(fileHandle, line);
    sizeNum = std::atoi(line.substr(23,line.length()).c_str());
    params.fileListMsk.resize(sizeNum);

    for (size_t i = 0; i < sizeNum; i++) {
        getline(fileHandle, line);
        params.fileListMsk[i] = line.substr(10,line.length());
    }

    fileHandle.close();
}

void msgTotalTime(sPickerParams& params, cTime& cTimer)
{
    std::stringstream   strStream;
    std::ofstream       fileHandle;

    strStream << "Total time:      " << std::fixed
              << std::setw(5) << params.mpiNumProcs << " nProc, "
              << std::setw(5) << params.nGPU << " nGPU, "
              << std::setw(5) << params.nCPU << " nCPU, "
              << "   " << cTimer.getElapsedTime() << std::endl;

#ifdef __GEM_SAVE_TIMING__
    fileHandle.open("time_gempicker_total.txt", std::ofstream::app |
                                                std::ofstream::out);
    fileHandle << strStream.str();
    fileHandle.close();
#endif

    std::cout << strStream.str();
}

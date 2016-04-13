/***********************************************************************
 *  File:       cStackRegTPS2.cpp
 *
 *  Purpose:    Implementation of a stack registration class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cStackRegTPS2.hpp"

#include "cTIFF.hpp"
#include "cRegTPS2.hpp"
#include "cRegTrackPoint2.hpp"
#include "cResultAligner.hpp"
#include "cFileSystem.hpp"

namespace gem {

template <typename T>
void cStackRegTPS2<T>::setParamsPointTracking(const cSize3& tplSizeHalf,
                                              const cSize3& searchRange,
                                              const cSize3& step,
                                              const cSize3& border)
{
    _tplSizeHalf = tplSizeHalf;
    _searchRange = searchRange;
           _step = step;
         _border = border;
}

template <typename T>
void cStackRegTPS2<T>::trackPointsUsingCorr(const std::string& dirDataIn,
                                            const std::string& dirDataTmp,
                                            const std::vector<std::string>& fileList,
                                            eCCmode mode)
{
    const std::string   funcName("void cStackRegTPS2<T>::trackPointsUsingCorr("
                                    "const std::string& dirDataIn, "
                                    "const std::string& dirDataTmp, "
                                    "const std::vector<std::string>& fileList, "
                                    "eCCmode mode)");

    require(fileList.size() > 2, funcName + ": not an image stack");

    cData3<T>               dataOriginal;
    cData3<T>               dataDistorted;

    cRegTrackPoint2<T>      objTrackPoint;
    std::string             dirSave;

    cTIFF<T>().read(dirDataIn + "/" + fileList[0], dataOriginal);

    objTrackPoint.genPointsRectangular(dataOriginal.getSize(),
                                       _border,
                                       _step);

    cResultAligner<T>().writePoints(objTrackPoint.getPtsOriginal(),
                                    dirDataTmp + "/ptsOriginal.txt");

    for (size_t i = 1; i < fileList.size(); i++) {

#ifndef NDEBUG
        std::cout << fileList[i] << std::endl;
#else
        std::cout << "\r     " << i+1 << "/" << num2str(fileList.size()) << std::flush;
#endif

        dirSave = dirDataTmp + "/" + cFileSystem(fileList[i]).getFileRoot();
        cFileSystem(dirSave).dirCreate();

        cTIFF<T>().read(dirDataIn + "/" + fileList[i], dataDistorted);

        require(dataOriginal.getSize() == dataDistorted.getSize(),
                funcName + ": images in stack do not have the same size");

        objTrackPoint.trackPointsUsingCorr(dataOriginal,
                                           dataDistorted,
                                           _tplSizeHalf,
                                           _searchRange,
                                           mode);

        cResultAligner<T>().writePoints  (objTrackPoint.getPtsDistorted(),
                                          dirSave + "/ptsDistorted.txt");
        cResultAligner<T>().writeCorrVals(objTrackPoint.getCorrValues(),
                                          dirSave + "/corrVals.txt");

        dataOriginal = dataDistorted;
    }
}

template <typename T>
void cStackRegTPS2<T>::computeMatrixW(const std::string& dirDataTmp,
                                      const std::vector<std::string>& fileList)
{
    const std::string   funcName("void cStackRegTPS2<T>::computeMatrixW("
                                    "const std::string& dirDataTmp, "
                                    "const std::vector<std::string>& fileList)");

    require(fileList.size() > 2, funcName + ": not an image stack");

    cRegTPS2<T>             objTPS2;
    std::string             dirSave;

    #pragma omp parallel for private(dirSave, objTPS2)
    for (size_t i = 1; i < fileList.size(); i++) {

#ifndef NDEBUG
        std::cout << fileList[i] << std::endl;
#endif

        dirSave = dirDataTmp + "/" + cFileSystem(fileList[i]).getFileRoot();

        objTPS2.readPtsOriginal (dirDataTmp  + "/ptsOriginal.txt");
        objTPS2.readPtsDistorted(dirSave + "/ptsDistorted.txt");

        objTPS2.computeMatrixW();

        objTPS2.writeTpsW(dirSave + "/tpsW.txt");
    }
}

template <typename T>
void cStackRegTPS2<T>::computeMappingAndInterp(const std::string& dirDataIn,
                                               const std::string& dirDataTmp,
                                               const std::string& dirDataOut,
                                               const std::vector<std::string>& fileList)
{
    const std::string   funcName("void cStackRegTPS2<T>::computeMappingAndInterp("
                                    "const std::string& dirDataIn, "
                                    "const std::string& dirDataTmp, "
                                    "const std::string& dirDataOut, "
                                    "const std::vector<std::string>& fileList)");

    require(fileList.size() > 2, funcName + ": not an image stack");

    cRegTPS2<T>             objTPS2;
    cData3<T>               dataOriginal, dataDistorted, dataCorrected;
    std::string             dirSave;

    cTIFF<T>().read(dirDataIn + "/" + fileList[0], dataOriginal);
    objTPS2.gen2Dgrid(dataOriginal.getSize());

#ifdef __GEM_USE_CUDA__
    cuRegTPS2<T>            objTPS2CUDA;
    cuData3<T>              dataDistortedCUDA, dataCorrectedCUDA;

    objTPS2CUDA.copyNDgridToGPU(objTPS2);
#endif

    for (size_t i = 1; i < fileList.size(); i++) {

#ifndef NDEBUG
        std::cout << fileList[i] << std::endl;
#else
        std::cout << "\r     " << i+1 << "/" << num2str(fileList.size()) << std::flush;
#endif

        dirSave = dirDataTmp + "/" + cFileSystem(fileList[i]).getFileRoot();

        cTIFF<T>().read(dirDataIn + "/" + fileList[i], dataDistorted);

        objTPS2.readPtsDistorted(dirSave + "/ptsDistorted.txt");
        objTPS2.readTpsW(dirSave + "/tpsW.txt");

#ifdef __GEM_USE_CUDA__
        objTPS2CUDA.copyDataToGPU(objTPS2);
        objTPS2CUDA.computeMapping();

        dataDistortedCUDA = dataDistorted;
        objTPS2CUDA.computeInterp(dataDistortedCUDA, dataCorrectedCUDA, INTER_LINEAR);
        dataCorrected = dataCorrectedCUDA;
#else
        objTPS2.computeMapping();

        objTPS2.computeInterp(dataDistorted, dataCorrected, INTER_LINEAR);
#endif

        // this requires memory
        //objTPS2.writeMaps(dirSave);

        cTIFF<T>().write(dataOriginal,  dirSave + "/dataOriginal.tif");
        cTIFF<T>().write(dataDistorted, dirSave + "/dataDistorted.tif");
        cTIFF<T>().write(dataCorrected, dirSave + "/dataCorrected.tif");

        //cTIFF().write(dataCorrected, dirDataOut + "/" + fileList[i]);

        dataOriginal = dataDistorted;
    }
}

// instantiation
template class cStackRegTPS2<float >;
template class cStackRegTPS2<double>;

} // namespace gem

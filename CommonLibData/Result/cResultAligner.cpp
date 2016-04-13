/***********************************************************************
 *  File:       cResultAligner.cpp
 *
 *  Purpose:    Implementation of gEMaligner's result structure
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cResultAligner.hpp"

#include "cFileSystem.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace gem {

template <typename T>
void cResultAligner<T>::swapPointsCoordinates(cData3<T>& points)
{
    const std::string   funcName("void cResultAligner<T>::swapPointsCoordinates(cData3<T>& points)");

    points.requireNonEmpty(funcName);

    require(points.getNrow() > 0 && points.getNcol() == 2 && points.getNsec() == 1,
            funcName + ": invalid points data");

    for (size_t i = 0; i < points.getNrow(); i++) {
        std::swap(points[2*i],points[2*i+1]);
    }
}

template <typename T>
void cResultAligner<T>::readCorrVals(const std::string& fileName,
                                     cData3<T>& corrs) const
{
    const std::string   funcName("void cResultAligner<T>::readCorrVals("
                                    "const std::string& fileName, "
                                    "cData3<T>& corrs) const");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".txt",
            funcName + ": invalid file name = " + fileName);

    std::ifstream           fileHandle;
    std::string             lineStr;
    size_t                  npoint, ind = 0;

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    std::getline(fileHandle, lineStr);

    std::stringstream       converterHeader(lineStr);
    converterHeader >> npoint;

    require(npoint != 0, funcName + ": no value detected in the file " + fileName);
    corrs.memReAlloc(cSize3(npoint,1,1));

    while (!fileHandle.eof()) {

        std::getline(fileHandle, lineStr);

        if (lineStr != "") {

            std::stringstream       converter(lineStr);

            converter >> corrs[ind++];
        }
    }

    fileHandle.close();
}

template <typename T>
void cResultAligner<T>::writeCorrVals(const cData3<T>& corrs,
                                      const std::string& fileName) const
{
    const std::string   funcName("void cResultAligner<T>::writeCorrVals("
                                    "const cData3<T>& corrs, "
                                    "const std::string& fileName) const");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".txt",
            funcName + ": invalid file name = " + fileName);

    unsigned int            precision = 4;
    std::ofstream           fileHandle;
    size_t                  npoint = corrs.getNrow();

    require(npoint > 0 && corrs.getNcol() == 1 && corrs.getNsec() == 1,
            funcName + ": invalid correlation data");

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    fileHandle << std::setprecision(precision) << std::fixed
               << std::setw( 8)           << npoint
               << std::setw( 8+precision) << "corr"
               << std::endl;

    for (size_t i = 0; i < npoint; i++) {

        fileHandle << std::setprecision(precision) << std::fixed
                   << std::setw(16+precision) << corrs[i]
                   << std::endl;
    }

    fileHandle.close();
}

template <typename T>
void cResultAligner<T>::readPoints(const std::string& fileName, cData3<T>& points) const
{
    const std::string   funcName("void cResultAligner<T>::readPoints("
                                    "const std::string& fileName, "
                                    "cData3<T>& points) const");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".txt",
            funcName + ": invalid file name = " + fileName);

    std::ifstream           fileHandle;
    std::string             lineStr;
    size_t                  npoint, ind = 0;
    T                       row, col;

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    std::getline(fileHandle, lineStr);

    std::stringstream       converterHeader(lineStr);
    converterHeader >> npoint;

    require(npoint != 0, funcName + ": no point detected in the file " + fileName);
    points.memReAlloc(cSize3(npoint,2,1));

    while (!fileHandle.eof()) {

        std::getline(fileHandle, lineStr);

        if (lineStr != "") {

            std::stringstream       converter(lineStr);

            converter >> row >> col;

            points[ind++] = row;
            points[ind++] = col;
        }
    }

    fileHandle.close();
}

template <typename T>
void cResultAligner<T>::writePoints(const cData3<T>& points, const std::string& fileName,
                                    unsigned int precision) const
{
    const std::string   funcName("void cResultAligner<T>::writePoints("
                                    "const cData3<T>& points, "
                                    "const std::string& fileName, "
                                    "unsigned int precision) const");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".txt",
            funcName + ": invalid file name = " + fileName);

    std::ofstream           fileHandle;
    size_t                  npoint = points.getNrow(), ind = 0;
    T                       row, col;

    require(npoint > 0 && points.getNcol() == 2 && points.getNsec() == 1,
            funcName + ": invalid points data");

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    fileHandle << std::setprecision(precision) << std::fixed
               << std::setw( 8)           << npoint
               << std::setw( 8+precision) << "row"
               << std::setw(16+precision) << "col"
               << std::endl;

    for (size_t i = 0; i < npoint; i++) {

        row = points[ind++];
        col = points[ind++];

        fileHandle << std::setprecision(precision) << std::fixed
                   << std::setw(16+precision) << row
                   << std::setw(16+precision) << col
                   << std::endl;
    }

    fileHandle.close();
}

template <typename T>
void cResultAligner<T>::readResultRigid(const std::string&          fileName,
                                        std::vector< cVector3<T> >& result) const
{
    const std::string   funcName("void cResultAligner<T>::readResultRigid("
                                    "const std::string& fileName, "
                                    "std::vector< cVector3<T> >& result) const");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".txt",
            funcName + ": invalid file name = " + fileName);

    std::ifstream           fileHandle;
    std::string             lineStr;
    size_t                  nRes, ind = 0;
    T                       row, col, corrVal;

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    std::getline(fileHandle, lineStr);

    std::stringstream       converterHeader(lineStr);
    converterHeader >> nRes;

    require(nRes != 0, funcName + ": no result detected in the file " + fileName);
    result.resize(nRes);

    while (!fileHandle.eof()) {

        std::getline(fileHandle, lineStr);

        if (lineStr != "") {

            std::stringstream       converter(lineStr);

            converter >> row >> col >> corrVal;

            result[ind][0] = row;
            result[ind][1] = col;
            result[ind][2] = corrVal;

            ind++;
        }
    }

    fileHandle.close();
}

template <typename T>
void cResultAligner<T>::writeResultRigid(const std::vector< cVector3<T> >& result,
                                         const std::string&                fileName,
                                         unsigned int                      precision) const
{
    const std::string   funcName("void cResultAligner<T>::writeResultRigid("
                                    "const std::vector< cVector3<T> >& result, "
                                    "const std::string& fileName, "
                                    "unsigned int precision) const");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".txt",
            funcName + ": invalid file name = " + fileName);

    std::ofstream           fileHandle;
    size_t                  nRes = result.size();

    require(nRes > 0, funcName + ": invalid result vector");

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    fileHandle << std::setprecision(precision) << std::fixed
               << std::setw( 8)           << nRes
               << std::setw( 8+precision) << "row"
               << std::setw(16+precision) << "col"
               << std::setw(16+precision) << "corrVal"
               << std::endl;

    for (size_t i = 0; i < nRes; i++) {

        fileHandle << std::setprecision(precision) << std::fixed
                   << std::setw(16+precision) << result[i][0]
                   << std::setw(16+precision) << result[i][1]
                   << std::setw(16+precision) << result[i][2]
                   << std::endl;
    }

    fileHandle.close();
}

template <typename T>
void cResultAligner<T>::checkResultRigid(const std::vector< cVector3<T> >& result,
                                         const cSize3&                     imSize,
                                         T                                 thrOffset,
                                         T                                 thrCorrVal) const
{
    const std::string   funcName("void cResultAligner<T>::checkResultRigid("
                                    "const std::vector< cVector3<T> >& result, "
                                    "const cSize3& imSize, "
                                    "T thrOffset, "
                                    "T thrCorrVal) const");

    size_t    nRes = result.size();
    std::vector<size_t>        vecOffset, vecCorrVal;

    require(nRes > 0,                          funcName + ": invalid result vector");
    require(imSize[0] > 0 && imSize[1] > 0,    funcName + ": invalid image size");
    require(thrOffset >= 0 && thrCorrVal >= 0 &&
            thrOffset <= 1 && thrCorrVal <= 1, funcName + ": invalid threshold values");

     vecOffset.clear();
    vecCorrVal.clear();
    for (size_t i = 0; i < nRes; i++) {

        if ((std::abs(result[i][0]) >= thrOffset * (T) imSize[0]) ||
            (std::abs(result[i][1]) >= thrOffset * (T) imSize[1])) {
            vecOffset.push_back(i);
        }

        if (result[i][2] <= thrCorrVal) {
            vecCorrVal.push_back(i);
        }
    }

    if (vecOffset.size()) {
        std::cout << "     large offset values > " << thrOffset << " x image size: ";

        for (size_t i = 0; i < vecOffset.size(); i++) {
            std::cout << vecOffset[i] << ", ";
        }
        std::cout << std::endl;
    }

    if (vecCorrVal.size()) {
        std::cout << "     small correlation values < " << thrCorrVal << ": ";

        for (size_t i = 0; i < vecCorrVal.size(); i++) {
            std::cout << vecCorrVal[i] << ", ";
        }
        std::cout << std::endl;
    }

    if (vecOffset.size() == 0 && vecCorrVal.size() == 0) {
        std::cout << "     no problem detected" << std::endl;
    }
}

// instantiation
template class cResultAligner<float >;
template class cResultAligner<double>;

} // namespace gem

/***********************************************************************
 *  File:       cResultPicker.cpp
 *
 *  Purpose:    Implementation of gEMpicker's result structure
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cResultPicker.hpp"

#include "macro.hpp"
#include "cFileSystem.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

namespace gem {

template <typename T>
std::ostream& operator<<(std::ostream& os, const sPickerResult<T>& result)
{
    std::cout << std::setprecision(4) << std::fixed
              << std::setw(4)  << result._numStr
              << std::setw(12) << result._corrVal
              << std::setw(10) << result._x
              << std::setw(10) << result._y
              << std::setw(10) << result._boxSize
              << std::setw(10) << result._tplIdx
              << std::setw(14) << result._angle2D
              << std::setw(30) << result._tplName
              << std::setw(30) << result._mskName
              << std::endl;

    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const cResultPicker<T>& result)
{
    std::cout << std::setprecision(4) << std::fixed
              << std::setw(4)  << "iPik"
              << std::setw(12) << "corrVal"
              << std::setw(10) << "X"
              << std::setw(10) << "Y"
              << std::setw(10) << "boxSize"
              << std::setw(10) << "iSch"
              << std::setw(14) << "angle2D"
              << std::setw(30) << "fileListSch[iSch]"
              << std::setw(30) << "fileListMsk[0/iSch]"
              << std::endl;

    size_t    nElem = result.size();
    for (size_t i = 0; i < nElem; i++) {

        std::cout << result.at(i);
    }

    return os;
}

template <typename T>
void cResultPicker<T>::read(const std::string& fileName)
{
    const std::string   funcName("void cResultPicker<T>::read("
                                    "const std::string& fileName)");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".txt",
            funcName + ": invalid file name = " + fileName);

    std::ifstream           fileHandle;
    std::string             lineStr;
    sPickerResult<T>        result;

    std::vector< sPickerResult<T> >::clear();

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    std::getline(fileHandle, lineStr);

    while (!fileHandle.eof()) {

        std::getline(fileHandle, lineStr);

        if (lineStr != "") {

            std::stringstream       converter(lineStr);

            converter >> result._numStr
                      >> result._corrVal
                      >> result._x
                      >> result._y
                      >> result._boxSize
                      >> result._tplIdx
                      >> result._angle2D
                      >> result._tplName
                      >> result._mskName;

            std::vector< sPickerResult<T> >::push_back(result);
        }
    }

    fileHandle.close();
}

template <typename T>
void cResultPicker<T>::write(const std::string& fileName) const
{
    const std::string   funcName("void cResultPicker<T>::write("
                                    "const std::string& fileName) const");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".txt",
            funcName + ": invalid file name = " + fileName);

    std::ofstream           fileHandle;
    sPickerResult<T>        result;
    size_t                  nElem = std::vector< sPickerResult<T> >::size();

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    fileHandle << std::setprecision(4) << std::fixed
               << std::setw(4)  << "iPik"
               << std::setw(12) << "corrVal"
               << std::setw(10) << "X"
               << std::setw(10) << "Y"
               << std::setw(10) << "boxSize"
               << std::setw(10) << "iSch"
               << std::setw(14) << "angle2D"
               << std::setw(30) << "fileListSch[iSch]"
               << std::setw(30) << "fileListMsk[0/iSch]"
               << std::endl;

    for (size_t i = 0; i < nElem; i++) {

        result = std::vector< sPickerResult<T> >::at(i);

        fileHandle << std::setprecision(4) << std::fixed
                   << std::setw(4)  << result._numStr
                   << std::setw(12) << result._corrVal
                   << std::setw(10) << result._x
                   << std::setw(10) << result._y
                   << std::setw(10) << result._boxSize
                   << std::setw(10) << result._tplIdx
                   << std::setw(14) << result._angle2D
                   << std::setw(30) << result._tplName
                   << std::setw(30) << result._mskName
                   << std::endl;
    }

    fileHandle.close();
}

template <typename T>
void cResultPicker<T>::writeBoxEMAN(const std::string& fileName, size_t refNrow) const
{
    const std::string   funcName("void cResultPicker<T>::writeBoxEMAN("
                                    "const std::string& fileName, "
                                    "size_t refNrow) const");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".box",
            funcName + ": invalid file name = " + fileName);

    std::ofstream           fileHandle;
    sPickerResult<T>        result;
    size_t                  nElem = std::vector< sPickerResult<T> >::size();

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    size_t boxSizeHalf;

    for (size_t i = 0; i < nElem; i++) {

        result      = std::vector< sPickerResult<T> >::at(i);
        boxSizeHalf = result._boxSize / 2 - 1;

        fileHandle << refNrow-1 - result._y - boxSizeHalf << "\t"
                   <<             result._x - boxSizeHalf << "\t"
                   << result._boxSize                     << "\t"
                   << result._boxSize
                   << std::endl;
    }

    fileHandle.close();
}

template <typename T>
size_t cResultPicker<T>::refineUsingDistance(size_t boxDist)
{
    const std::string   funcName("size_t cResultPicker<T>::refineUsingDistance(size_t boxDist)");

    require(std::vector< sPickerResult<T> >::size() > 0, funcName + ": empty cResultPicker object");

    std::vector<size_t>     ind;
    sPickerResult<T>        ri, rj;
    double                  distij;
    size_t                  nElemBefore = std::vector< sPickerResult<T> >::size();
    size_t                  nElemAfter;
    cString                 numStr;

    // determine too-close particles and store their identity in a vector
    for (size_t i = 0; i < nElemBefore-1; i++) {
        for (size_t j = i+1; j < nElemBefore; j++) {

            ri = std::vector< sPickerResult<T> >::at(i);
            rj = std::vector< sPickerResult<T> >::at(j);

            distij = dist( (double) ri._x, (double) ri._y,
                           (double) rj._x, (double) rj._y );

            if (distij < (double) boxDist) {
                ind.push_back(i+1);
                ind.push_back(j+1);
            }
        }
    }

    // remove redundant identities from vector
    std::sort(ind.begin(), ind.end());
    ind.erase(std::unique(ind.begin(), ind.end() ), ind.end());

    // remove particles using remaining identities
    for (size_t i = ind.size(); i > 0; i--) {
        std::vector< sPickerResult<T> >::erase(std::vector< sPickerResult<T> >::begin()+ind.at(i-1)-1);
    }

    // correct index
    nElemAfter = std::vector< sPickerResult<T> >::size();
    for (size_t i = 0; i < nElemAfter; i++) {
        numStr = num2str(i+1);
        numStr.insertBegin('0',4);

        std::vector< sPickerResult<T> >::at(i)._numStr = numStr;
    }

    // return the number of removed particles
    return nElemBefore - nElemAfter;
}

template <typename T>
size_t cResultPicker<T>::refineUsingHighThreshold(float thresholdHigh)
{
    const std::string   funcName("size_t cResultPicker<T>::refineUsingHighThreshold(float thresholdHigh)");

    require(std::vector< sPickerResult<T> >::size() > 0, funcName + ": empty cResultPicker object");

    size_t      removed = 0;
    size_t      nElem = std::vector< sPickerResult<T> >::size();

    // remove particles with correlation value over the threshold
    for (size_t i = 0; i < nElem; i++) {
        if (std::abs(std::vector< sPickerResult<T> >::at(i)._corrVal) > thresholdHigh) {
            std::vector< sPickerResult<T> >::erase(std::vector< sPickerResult<T> >::begin()+i);
            removed++;
            i--;
        }
    }

    // return the number of removed particles
    return removed;
}

// instantiation
template std::ostream& operator<<(std::ostream& os, const sPickerResult<float >& result);
template std::ostream& operator<<(std::ostream& os, const sPickerResult<double>& result);

template std::ostream& operator<<(std::ostream& os, const cResultPicker<float >& result);
template std::ostream& operator<<(std::ostream& os, const cResultPicker<double>& result);

template class cResultPicker<float >;
template class cResultPicker<double>;

} // namespace gem

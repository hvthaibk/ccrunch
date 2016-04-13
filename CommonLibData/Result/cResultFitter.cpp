/***********************************************************************
 *  File:       cResultFitter.cpp
 *
 *  Purpose:    Implementation of gEMfitter's result structure
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cResultFitter.hpp"

#include "macro.hpp"
#include "cFileSystem.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace gem {

template <typename T>
std::ostream& operator<<(std::ostream& os, const sFitterResult<T>& result)
{
    std::cout << std::setprecision(4) << std::fixed
              << std::setw(6)  << result._numStr
              << std::setw(18) << result._corrVal
              << std::setw(12) << result._pos[0]
              << std::setw(12) << result._pos[1]
              << std::setw(12) << result._pos[2]
              << std::setw(12) << result._angle[0]
              << std::setw(12) << result._angle[1]
              << std::setw(12) << result._angle[2]
              << std::endl;

    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const cResultFitter<T>& result)
{
    std::cout << std::setprecision(4) << std::fixed
              << std::setw(6)  << "numStr"
              << std::setw(18) << "corrVal"
              << std::setw(12) << "deltaRow"
              << std::setw(12) << "deltaCol"
              << std::setw(12) << "deltaSec"
              << std::setw(12) << "Alpha"
              << std::setw(12) << "Beta"
              << std::setw(12) << "Gamma"
              << std::endl;

    size_t    nElem = result.size();
    for (size_t i = 0; i < nElem; i++) {

        std::cout << result.at(i);
    }

    return os;
}

template <typename T>
void cResultFitter<T>::read(const std::string& fileName)
{
    const std::string   funcName("void cResultFitter<T>::read("
                                    "const std::string& fileName)");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".txt",
            funcName + ": invalid file name = " + fileName);

    std::ifstream           fileHandle;
    std::string             lineStr;
    sFitterResult<T>        result;

    std::vector< sFitterResult<T> >::clear();

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    std::getline(fileHandle, lineStr);

    while (!fileHandle.eof()) {

        std::getline(fileHandle, lineStr);

        if (lineStr != "") {

            std::stringstream       converter(lineStr);

            converter >> result._numStr
                      >> result._corrVal
                      >> result._pos[0]
                      >> result._pos[1]
                      >> result._pos[2]
                      >> result._angle[0]
                      >> result._angle[1]
                      >> result._angle[2];

            std::vector< sFitterResult<T> >::push_back(result);
        }
    }

    fileHandle.close();
}

template <typename T>
void cResultFitter<T>::write(const std::string& fileName) const
{
    const std::string   funcName("void cResultFitter<T>::write("
                                    "const std::string& fileName) const");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".txt",
            funcName + ": invalid file name = " + fileName);

    std::ofstream           fileHandle;
    sFitterResult<T>        result;
    size_t                  nElem = std::vector< sFitterResult<T> >::size();

    warning(nElem > 0, funcName + ": no result to write");

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    fileHandle << std::setprecision(4) << std::fixed
               << std::setw(6)  << "numStr"
               << std::setw(18) << "corrVal"
               << std::setw(12) << "deltaRow"
               << std::setw(12) << "deltaCol"
               << std::setw(12) << "deltaSec"
               << std::setw(12) << "Alpha"
               << std::setw(12) << "Beta"
               << std::setw(12) << "Gamma"
               << std::endl;

    for (size_t i = 0; i < nElem; i++) {

        result = std::vector< sFitterResult<T> >::at(i);

        fileHandle << std::setprecision(4) << std::fixed
                   << std::setw(6)  << result._numStr
                   << std::setw(18) << result._corrVal
                   << std::setw(12) << result._pos[0]
                   << std::setw(12) << result._pos[1]
                   << std::setw(12) << result._pos[2]
                   << std::setw(12) << result._angle[0]
                   << std::setw(12) << result._angle[1]
                   << std::setw(12) << result._angle[2]
                   << std::endl;
    }

    fileHandle.close();
}

// instantiation
template std::ostream& operator<<(std::ostream& os, const sFitterResult<float >& result);
template std::ostream& operator<<(std::ostream& os, const sFitterResult<double>& result);

template std::ostream& operator<<(std::ostream& os, const cResultFitter<float >& result);
template std::ostream& operator<<(std::ostream& os, const cResultFitter<double>& result);

template class cResultFitter<float >;
template class cResultFitter<double>;

} // namespace gem

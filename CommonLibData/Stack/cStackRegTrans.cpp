/***********************************************************************
 *  File:       cStackRegTrans.cpp
 *
 *  Purpose:    Implementation of a stack registration class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cStackRegTrans.hpp"

#include "cSystem.hpp"
#include "cTIFF.hpp"
#include "cRegCorr.hpp"
#include "cTransSimilarity.hpp"

namespace gem {

template <typename T>
void cStackRegTrans<T>::memFree(void)
{
    _corrRes.clear();
}

template <typename T>
std::vector< cVector3<T> >& cStackRegTrans<T>::getResult(void)
{
    const std::string   funcName("std::vector< cVector3<T> >& cStackRegTrans<T>::getResult(void)");

    require(_corrRes.size() > 0, funcName + ": correlation result vector is empty");

    return _corrRes;
}

template <typename T>
void cStackRegTrans<T>::setResult(const std::vector< cVector3<T> >& corrRes)
{
    const std::string   funcName("void cStackRegTrans<T>::setResult("
                                    "const std::vector< cVector3<T> >& corrRes)");

    require(corrRes.size() > 0, funcName + ": invalid result vector");

    _corrRes= corrRes;
}

template <typename T>
void cStackRegTrans<T>::setROIs(const cSize3& tplSize, const cSize3& tplOffset,
                                const cSize3& refSize, const cSize3& refOffset)
{
    const std::string   funcName("void cStackRegTrans<T>::setROIs("
                                    "const cSize3& tplSize, "
                                    "const cSize3& tplOffset, "
                                    "const cSize3& refSize, "
                                    "const cSize3& refOffset)");

    _tplSize   = tplSize;
    _tplOffset = tplOffset;

    _refSize   = refSize;
    _refOffset = refOffset;

    require(_tplOffset >= _refOffset, funcName + ": invalid offset values");
    require(_tplOffset+_tplSize <= _refOffset+_refSize,
            funcName + ": invalid offset+size values");
}

template <typename T>
void cStackRegTrans<T>::computeOffsetsUsingCorr(const std::string& dirDataIn,
        const std::vector<std::string>& fileList, eCCmode mode)
{
    const std::string   funcName("void cStackRegTrans<T>::computeOffsetsUsingCorr("
                                    "const std::string& dirDataIn, "
                                    "const std::vector<std::string>& fileList, "
                                    "eCCmode mode)");

    require(fileList.size() > 2, funcName + ": not an image stack");

    cTIFF<T>            objTIFF;
    cRegCorr<T>         objCorr;
    cData3XCorr<T>      dataBefore, dataAfter;
    cData3XCorr<T>      dataTpl, dataRef, dataRes;
    cSize3              resDelta(0,0,0);

    objTIFF.read(dirDataIn + "/" + fileList[0], dataBefore);

    _corrRes.resize(fileList.size());
    _corrRes[0] = cVector3<T>(0,0,1);   // first section

    for (size_t i = 1; i < fileList.size(); i++) {

#ifndef NDEBUG
        std::cout << fileList[i] << std::endl;
#else
        std::cout << "\r     " << i+1 << "/" << num2str(fileList.size()) << std::flush;
#endif

        objTIFF.read(dirDataIn + "/" + fileList[i], dataAfter);

        require(dataBefore.getSize() == dataAfter.getSize(),
                funcName + ": images in stack do not have the same size");

        dataTpl.opCrop(dataBefore, _tplSize, _tplOffset);
        dataRef.opCrop(dataAfter,  _refSize, _refOffset);

        dataBefore = dataAfter;

        switch (mode) {
            case CC_MODE_XCORR:
                objCorr.computeXCC(dataTpl, dataRef, dataRes, XCORR_RES_VALID);
                break;
            case CC_MODE_ECORR:
                objCorr.computeECC(dataTpl, dataRef, dataRes, XCORR_RES_VALID, false);
                break;
            case CC_MODE_NCORR:
                objCorr.computeNCC(dataTpl, dataRef, dataRes, XCORR_RES_VALID, false);

                require(dataRes.isFinite(), funcName + ": invalid correlation results");
                break;
            default:
                ERROR(funcName, "unsupported correlation mode");
        }

        resDelta.ind2sub(dataRes.getSize(), dataRes.getMaxIndx());

        // store offset and correlation values
        _corrRes[i][0] = (T) resDelta[0] - (T) _tplOffset[0] + (T) _refOffset[0];
        _corrRes[i][1] = (T) resDelta[1] - (T) _tplOffset[1] + (T) _refOffset[1];
        _corrRes[i][2] = dataRes.getMax();
    }
}

template <typename T>
void cStackRegTrans<T>::alignSections(const std::string& dirDataIn,
                                      const std::string& dirDataOut,
                                      const std::vector<std::string>& fileList,
                                      bool  extend)
{
    const std::string   funcName("void cStackRegTrans<T>::alignSections("
                                    "const std::string& dirDataIn, "
                                    "const std::string& dirDataOut, "
                                    "const std::vector<std::string>& fileList, "
                                    "bool extend)");

    require(fileList.size() > 2, funcName + ": not an image stack");
    require(_corrRes.size() == fileList.size(),
            funcName + ": distances need to be computed a priori");

    cTIFF<T>                objTIFF;
    cTransSimilarity<T>     objTrans;
    cData3<T>               data, dataExt;
    cData3<T>               offRow, offCol;
    T                       offRowMin(0), offRowMax(0);
    T                       offColMin(0), offColMax(0);

    offRow.memReAlloc(cSize3(_corrRes.size(),1,1));
    offCol.memReAlloc(cSize3(_corrRes.size(),1,1));
    offRow[0] = 0;
    offCol[0] = 0;

    for (size_t i = 1; i < fileList.size(); i++) {

        offRow[i] = offRow[i-1] + _corrRes[i][0];
        offCol[i] = offCol[i-1] + _corrRes[i][1];
    }

    offRowMin = offRow.getMin();
    offRowMax = offRow.getMax();
    offColMin = offCol.getMin();
    offColMax = offCol.getMax();

    objTIFF.read (dirDataIn + "/" + fileList[0], data);
    if (extend) {
        dataExt.memReAllocZero(cSize3((size_t) ((T) data.getNrow()-offRowMin+offRowMax),
                                      (size_t) ((T) data.getNcol()-offColMin+offColMax),
                                      1));
        dataExt.opReplace(data, cSize3((size_t) offRowMax,
                                       (size_t) offColMax,
                                       0));
        objTIFF.write(dataExt, dirDataOut + "/" + fileList[0]);
    }
    else {
        objTIFF.write(data, dirDataOut + "/" + fileList[0]);
    }

    for (size_t i = 1; i < fileList.size(); i++) {

#ifndef NDEBUG
        std::cout << fileList[i] << std::endl;
#else
        std::cout << "\r     " << i+1 << "/" << num2str(fileList.size()) << std::flush;
#endif

        objTIFF.read(dirDataIn + "/" + fileList[i], data);

        if (extend) {
            dataExt.memSetZero();
            dataExt.opReplace(data, cSize3((size_t)(offRowMax-offRow[i]),
                                           (size_t)(offColMax-offCol[i]),
                                           0));

            objTIFF.write(dataExt, dirDataOut + "/" + fileList[i]);
        }
        else {
            objTrans.translate(data, cVector3<T>(-offRow[i],-offCol[i],0));

            objTIFF.write(data, dirDataOut + "/" + fileList[i]);
        }
    }
}

// instantiation
template class cStackRegTrans<float >;
template class cStackRegTrans<double>;

} // namespace gem

/***********************************************************************
 *  File:       cRegTPS2.cpp
 *
 *  Purpose:    Implementation of a TPS2 registration class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cRegTPS2.hpp"

#include "macro.hpp"
#include "cNDgrid.hpp"
#include "cResultAligner.hpp"
#include "cMRC.hpp"

#include "cEigen.hpp"
#include "cArmadillo.hpp"

#include <cmath>

namespace gem {

template <typename T>
void cRegTPS2<T>::memFree(void)
{
            _tpsL.memFree();
            _tpsW.memFree();
     _ptsOriginal.memFree();
    _ptsDistorted.memFree();
               _X.memFree();
               _Y.memFree();
            _Xmap.memFree();
            _Ymap.memFree();
}

template <typename T>
void cRegTPS2<T>::gen2Dgrid(cSize3 size)
{
    const std::string   funcName("void cRegTPS2<T>::gen2Dgrid(cSize3 size)");

    require(size[0] > 1 && size[1] > 1 && size[2] == 1,
            funcName + ": size is not 2D");

    cData3<T>   xVec, yVec;

    xVec.linspace(0, (T) size[1]-1, size[1]);
    yVec.linspace(0, (T) size[0]-1, size[0]);

    gen2Dgrid(xVec, yVec);
}

template <typename T>
void cRegTPS2<T>::gen2Dgrid(cData3<T>& xVec, cData3<T>& yVec)
{
    const std::string   funcName("void cRegTPS2<T>::gen2Dgrid("
                                    "cData3<T>& xVec, cData3<T>& yVec)");

    cNDgrid<T>().genXY(xVec, yVec, _X, _Y);
}

template <typename T>
void cRegTPS2<T>::computeMatrixL(void)
{
    const std::string   funcName("void cRegTPS2<T>::computeMatrixL(void)");

    size_t      npoint = _ptsDistorted.getNrow();
    T           valTmp;
    cData3<T>   tpsK, tpsP;
    cData3<T>   Xd, Yd;

    Xd.getSlideCol(_ptsDistorted,0);
    Yd.getSlideCol(_ptsDistorted,1);

    // tpsK
    tpsK.memReAlloc(cSize3(npoint, npoint, 1));
    for (size_t i = 0; i < npoint; i++) {
        for (size_t j = 0; j < npoint; j++) {

            valTmp = distSqr(Xd[i],Yd[i],Xd[j],Yd[j]);

            if (valTmp != 0) { tpsK[sub2ind(i,j,npoint)] = valTmp * std::log(valTmp); }
            else             { tpsK[sub2ind(i,j,npoint)] = 0;                         }
        }
    }

    // tpsP
    tpsP.memReAlloc(cSize3(npoint,3,1));
    tpsP.memSetVal(1);
    tpsP.opReplace(_ptsDistorted,cSize3(0,1,0));

    // tpsL
    _tpsL.memReAllocZero(cSize3(npoint+3,npoint+3,1));
    _tpsL.opReplace(tpsK,cSize3(0,0,0));
    _tpsL.opReplace(tpsP,cSize3(0,npoint,0));
     tpsP.opPermute(PERMUTE2D);
    _tpsL.opReplace(tpsP,cSize3(npoint,0,0));
}

template <typename T>
void cRegTPS2<T>::computeMatrixW(void)
{
    const std::string   funcName("void cRegTPS2<T>::computeMatrixW(void)");

    _ptsOriginal.requireNonEmpty(funcName + ": _ptsOriginal cannot be empty");
    require(_ptsOriginal.getNrow() > 0 && _ptsOriginal.getNcol() == 2,
            funcName + ": _ptsOriginal data is incorrect");

    _ptsDistorted.requireNonEmpty(funcName + ": _ptsDistorted cannot be empty");
    require(_ptsDistorted.getNrow() > 0 && _ptsDistorted.getNcol() == 2,
            funcName + ": _ptsDistorted data is incorrect");

    require(_ptsOriginal.getSize() == _ptsDistorted.getSize(),
            funcName + ": _ptsOriginal and _ptsDistorted are not compatible");

    size_t      npoint = _ptsDistorted.getNrow();
    cData3<T>   tpsY;

    computeMatrixL();

    tpsY.memReAllocZero(cSize3(npoint+3,2,1));
    tpsY.opReplace(_ptsOriginal,cSize3(0,0,0));

    // use Armadillo/Eigen to solve Ax=b
    //cArmadillo<T>().solveAxb(_tpsL,tpsY,_tpsW);
    cEigen<T>().solveAxb(_tpsL,tpsY,_tpsW);
}

template <typename T>
void cRegTPS2<T>::computeMapping(void)
{
    const std::string   funcName("void cRegTPS2<T>::computeMapping(void)");

            _tpsW.requireNonEmpty(funcName + ": _tpsW cannot be empty");
               _X.requireNonEmpty(funcName + ": _X cannot be empty");
               _Y.requireNonEmpty(funcName + ": _Y cannot be empty");
    _ptsDistorted.requireNonEmpty(funcName + ": _ptsDistorted cannot be empty");

    require(_X.getSize() == _Y.getSize(),
            funcName + ": _X, _Y must have the same size");
    require(_tpsW.getNcol() == 2 && _tpsW.getNsec() == 1,
            funcName + ": _tpsW has incorrect size");
    require(_tpsW.getSize() == _ptsDistorted.getSize() + cSize3(3,0,0),
            funcName + ": _tpsW and _ptsDistorted are incompatible");

    size_t      npoint = _ptsDistorted.getNrow();
    size_t      npixel = _X.getNelement();
    size_t      i, j;
    T           valTmp, tpsK;
    cData3<T>   Xo, Yo, Xd, Yd;

    Xo.getSlideCol(_tpsW,0);
    Yo.getSlideCol(_tpsW,1);
    Xd.getSlideCol(_ptsDistorted,0);
    Yd.getSlideCol(_ptsDistorted,1);

    _Xmap.memReAlloc(_X.getSize());
    _Ymap.memReAlloc(_X.getSize());

    _Xmap.memSetVal(Xo[npoint]);
    _Xmap +=   _X * Xo[npoint+1];
    _Xmap +=   _Y * Xo[npoint+2];

    _Ymap.memSetVal(Yo[npoint]);
    _Ymap +=   _X * Yo[npoint+1];
    _Ymap +=   _Y * Yo[npoint+2];

    #pragma omp parallel for private(i,valTmp,tpsK)
    for (j = 0; j < npixel; j++) {
        for (i = 0; i < npoint; i++) {

            valTmp = distSqr(Xd[i],Yd[i],_X[j],_Y[j]);

            if (valTmp != 0) { tpsK = valTmp * std::log(valTmp); }
            else             { tpsK = 0;                         }

            _Xmap[j] += tpsK * Xo[i];
            _Ymap[j] += tpsK * Yo[i];
        }
    }
}

template <typename T>
void cRegTPS2<T>::computeInterp(const cData3<T>& dataDistorted, cData3<T>& dataCorrected, eInter inter)
{
    const std::string   funcName("void cRegTPS2<T>::computeInterp("
                                    "const cData3<T>& dataDistorted, "
                                    "cData3<T>& dataCorrected, eInter inter)");

    cTransSimilarity<T>().interp(_Ymap, _Xmap, dataDistorted, dataCorrected, inter);
}

template <typename T>
void cRegTPS2<T>::readPtsOriginal(const std::string& fileName)
{
    const std::string   funcName("void cRegTPS2<T>::readPtsOriginal("
                                    "const std::string& fileName)");

    cResultAligner<T>().readPoints(fileName,  _ptsOriginal);
    cResultAligner<T>().swapPointsCoordinates(_ptsOriginal);

    _ptsOriginal.requireNonEmpty(funcName);
}

template <typename T>
void cRegTPS2<T>::readPtsDistorted(const std::string& fileName)
{
    const std::string   funcName("void cRegTPS2<T>::readPtsDistorted("
                                    "const std::string& fileName)");

    cResultAligner<T>().readPoints(fileName,  _ptsDistorted);
    cResultAligner<T>().swapPointsCoordinates(_ptsDistorted);

    _ptsDistorted.requireNonEmpty(funcName);

}

template <typename T>
void cRegTPS2<T>::readTpsW(const std::string& fileName)
{
    const std::string   funcName("void cRegTPS2<T>::readTpsW("
                                    "const std::string& fileName)");

    cResultAligner<T>().readPoints(fileName, _tpsW);

    _tpsW.requireNonEmpty(funcName);
}

template <typename T>
void cRegTPS2<T>::writeTpsW(const std::string& fileName)
{
    const std::string   funcName("void cRegTPS2<T>::writeTpsW("
                                    "const std::string& fileName)");

    _tpsW.requireNonEmpty(funcName);

    cResultAligner<T>().writePoints(_tpsW, fileName, 16);
}

template <typename T>
void cRegTPS2<T>::readMaps(const std::string& dirName)
{
    const std::string   funcName("void cRegTPS2<T>::readMaps("
                                    "const std::string& dirName)");

    cMRC<T>().read(dirName + "/Xmap.mrc", _Xmap);
    cMRC<T>().read(dirName + "/Ymap.mrc", _Ymap);

    _Xmap.requireNonEmpty(funcName);
    _Ymap.requireNonEmpty(funcName);
}

template <typename T>
void cRegTPS2<T>::writeMaps(const std::string& dirName)
{
    const std::string   funcName("void cRegTPS2<T>::writeMaps("
                                    "const std::string& dirName)");

    _Xmap.requireNonEmpty(funcName);
    _Ymap.requireNonEmpty(funcName);

    cMRC<T>().write(_Xmap, dirName + "/Xmap.mrc");
    cMRC<T>().write(_Ymap, dirName + "/Ymap.mrc");
}

// instantiation
template class cRegTPS2<float >;
template class cRegTPS2<double>;

} // namespace gem

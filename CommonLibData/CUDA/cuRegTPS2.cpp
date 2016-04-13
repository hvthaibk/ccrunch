/***********************************************************************
 *  File:       cuRegTPS2.cpp
 *
 *  Purpose:    Implementation of a TPS2 registration class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuRegTPS2.hpp"

#include "registration.cuh"
#include "cuTransSimilarity.hpp"

namespace gem {

template <typename T>
void cuRegTPS2<T>::memFree(void)
{
            _tpsW.memFree();
    _ptsDistorted.memFree();
               _X.memFree();
               _Y.memFree();
            _Xmap.memFree();
            _Ymap.memFree();
}

template <typename T>
void cuRegTPS2<T>::copyNDgridToGPU(const cRegTPS2<T>& other)
{
    _X = other._X;
    _Y = other._Y;
}

template <typename T>
void cuRegTPS2<T>::copyDataToGPU(const cRegTPS2<T>& other)
{
            _tpsW = other._tpsW;
    _ptsDistorted = other._ptsDistorted;
}

template <typename T>
void cuRegTPS2<T>::copyDataFromGPU(cRegTPS2<T>& other) const
{
    other._Xmap = _Xmap;
    other._Ymap = _Ymap;
}

template <typename T>
void cuRegTPS2<T>::computeMapping(void)
{
    const std::string   funcName("void cuRegTPS2<T>::computeMapping(void)");

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
    cData3<T>    Xo,  Yo,  Xd,  Yd;
    cuData3<T>  _Xo, _Yo, _Xd, _Yd;

    Xo.getSlideCol(_tpsW,0);
    Yo.getSlideCol(_tpsW,1);
    Xd.getSlideCol(_ptsDistorted,0);
    Yd.getSlideCol(_ptsDistorted,1);

    _Xo = Xo;    _Yo = Yo;
    _Xd = Xd;    _Yd = Yd;

    _Xmap.memReAlloc(_X.getSize());
    _Ymap.memReAlloc(_X.getSize());

    cuda_regTPS2_compute_mapping(  _Xo.getAddrData(),   _Yo.getAddrData(),
                                   _Xd.getAddrData(),   _Yd.getAddrData(),
                                    _X.getAddrData(),    _Y.getAddrData(),
                                 _Xmap.getAddrData(), _Ymap.getAddrData(),
                                 npoint, npixel);
}

template <typename T>
void cuRegTPS2<T>::computeInterp(const cuData3<T>& dataDistorted, cuData3<T>& dataCorrected, eInter inter)
{
    const std::string   funcName("void cuRegTPS2<T>::computeInterp("
                                    "const cuData3<T>& dataDistorted, "
                                    "cuData3<T>& dataCorrected, eInter inter)");

    cuTransSimilarity<T>().interp(_Ymap, _Xmap, dataDistorted, dataCorrected, inter);
}

// instantiation
template class cuRegTPS2<float >;
template class cuRegTPS2<double>;

} // namespace gem

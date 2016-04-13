/***********************************************************************
 *  File:       cRegTPS2.hpp
 *
 *  Purpose:    Header file for a TPS2 registration class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CREG_TPS2_HPP__
#define __GEM_CREG_TPS2_HPP__

#include "cData3.hpp"
#include "cTransSimilarity.hpp"

#ifdef __GEM_USE_CUDA__
#include "cuRegTPS2.hpp"
#endif

namespace gem {

template <typename T>
class cRegTPS2
{
#ifdef __GEM_USE_CUDA__
template<typename U> friend class cuRegTPS2;
#endif

private:
    cData3<T>   _tpsL, _tpsW;
    cData3<T>   _ptsOriginal, _ptsDistorted;
    cData3<T>   _X, _Y, _Xmap, _Ymap;

public:
     cRegTPS2(void) {};
    ~cRegTPS2()     {};

    void    memFree(void);

    // set and get internal data
    void    setTpsL        (const cData3<T>& tpsL        ) {         _tpsL = tpsL;         };
    void    setTpsW        (const cData3<T>& tpsW        ) {         _tpsW = tpsW;         };
    void    setX           (const cData3<T>& X           ) {            _X = X;            };
    void    setY           (const cData3<T>& Y           ) {            _Y = Y;            };
    void    setPtsOriginal (const cData3<T>& ptsOriginal ) { _ptsOriginal  = ptsOriginal;  };
    void    setPtsDistorted(const cData3<T>& ptsDistorted) { _ptsDistorted = ptsDistorted; };

    void    getTpsL(cData3<T>& tpsL) { tpsL = _tpsL; };
    void    getTpsW(cData3<T>& tpsW) { tpsW = _tpsW; };
    void    getXmap(cData3<T>& Xmap) { Xmap = _Xmap; };
    void    getYmap(cData3<T>& Ymap) { Ymap = _Ymap; };

    // 2D grids
    void    gen2Dgrid(cSize3 size);
    void    gen2Dgrid(cData3<T>& xVec, cData3<T>& yVec);

    // compute TPS
    void    computeMatrixL(void);
    void    computeMatrixW(void);
    void    computeMapping(void);
    void    computeInterp (const cData3<T>& dataDistorted,
                                 cData3<T>& dataCorrected,
                           eInter inter = INTER_LINEAR);

    // save/load computed data
    void    readPtsOriginal (const std::string& fileName);
    void    readPtsDistorted(const std::string& fileName);

    void    readTpsW (const std::string& fileName);
    void    writeTpsW(const std::string& fileName);

    void    readMaps (const std::string&  dirName);
    void    writeMaps(const std::string&  dirName);
};

typedef     cRegTPS2<float >    cRegTPS2f;
typedef     cRegTPS2<double>    cRegTPS2d;

} // namespace gem

#endif

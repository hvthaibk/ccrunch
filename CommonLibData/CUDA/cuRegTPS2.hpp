/***********************************************************************
 *  File:       cuRegTPS2.hpp
 *
 *  Purpose:    Header file for a TPS2 registration class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUREG_TPS2_HPP__
#define __GEM_CUREG_TPS2_HPP__

#include "cRegTPS2.hpp"
#include "cuData3.hpp"

namespace gem {

template <typename T> class cRegTPS2;

template <typename T>
class cuRegTPS2
{
private:
    cData3<T>   _tpsW;
    cData3<T>   _ptsDistorted;
    cuData3<T>  _X, _Y, _Xmap, _Ymap;

public:
     cuRegTPS2(void) {};
    ~cuRegTPS2()     {};

    void    memFree(void);

    // import and export data from/to GPU
    void    copyNDgridToGPU(const cRegTPS2<T>& other);
    void    copyDataToGPU  (const cRegTPS2<T>& other);
    void    copyDataFromGPU(cRegTPS2<T>& other) const;

    // compute TPS
    void    computeMapping(void);
    void    computeInterp (const cuData3<T>& dataDistorted,
                                 cuData3<T>& dataCorrected,
                           eInter inter = INTER_LINEAR);
};

typedef     cuRegTPS2<float >   cuRegTPS2f;
typedef     cuRegTPS2<double>   cuRegTPS2d;

} // namespace gem

#endif

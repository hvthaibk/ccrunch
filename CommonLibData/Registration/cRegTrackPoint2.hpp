/***********************************************************************
 *  File:       cRegTrackPoint2.hpp
 *
 *  Purpose:    Header file for a 2D point tracking class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CREG_TRACKPOINT2_HPP__
#define __GEM_CREG_TRACKPOINT2_HPP__

#include "xcorr.hpp"
#include "cData3.hpp"

namespace gem {

template <typename T>
class cRegTrackPoint2
{
private:
    cData3<T>   _ptsOriginal;
    cData3<T>   _ptsDistorted;
    cData3<T>   _corrValue;

public:
     cRegTrackPoint2(void) {};
    ~cRegTrackPoint2()     {};

    void    memFree(void);

    cData3<T>&  getPtsOriginal (void);
    cData3<T>&  getPtsDistorted(void);
    cData3<T>&  getCorrValues  (void);

    void        setPointsOriginal(const cData3<T>& points);

    // points generation
    void    genPointsRectangular(const cSize3& size,
                                 const cSize3& border,
                                 const cSize3& step);

    // points tracking
    void    trackPointsUsingCorr(const cData3<T>& dataOriginal,
                                 const cData3<T>& dataDistorted,
                                 const cSize3&    tplSizeHalf,
                                 const cSize3&    searchRange,
                                 eCCmode mode = CC_MODE_NCORR);
};

} // namespace gem

#endif

/***********************************************************************
 *  File:       cTransDistance.hpp
 *
 *  Purpose:    Header file for a distance transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CTRANS_DISTANCE_HPP__
#define __GEM_CTRANS_DISTANCE_HPP__

#include "cData3.hpp"

namespace gem {

template <typename T>
class cTransDistance
{
private:

public:
     cTransDistance(void) {};
    ~cTransDistance()     {};

    // distance
    void    distanceGeneric(const cData3<T>& dataSrc, cData3<T>& dataDst);
    void    distanceNormal (const cData3<T>& dataSrc, cData3<T>& dataDst);

    // Voronoi
    void    voronoi(const cData3<T>& dataSrc, const cArray3<size_t>& labelSrc,
                          cData3<T>& dataDst,       cArray3<size_t>& labelDst);

    // watershed
    void    watershed(const cData3<T>& dataSrc, cArray3<size_t>& labelDst,
                      T cutoff, T seedDistMin,
                      eWaterDir water = WATER_DIR_RISE);

    void    watershedReverseLabel(cArray3<size_t>& labelSrc);

    // label
    void    label(const cData3<T>& dataSrc, cArray3<size_t>& labelDst);
};

} // namespace gem

#endif

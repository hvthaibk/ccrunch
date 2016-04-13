/***********************************************************************
 *  File:       cTransThreshold.hpp
 *
 *  Purpose:    Header file for a thresholding transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CTRANS_THRESHOLD_HPP__
#define __GEM_CTRANS_THRESHOLD_HPP__

#include "cData3.hpp"

namespace gem {

template <typename T>
class cTransThreshold
{
private:

public:
     cTransThreshold(void) {};
    ~cTransThreshold()     {};

    // thresholding
    T       getCutOffValue(const cData3<T>& dataSrc, T voxelVolume, T mapVolume,
                           eWaterDir water = WATER_DIR_FALL);
    void    setValueBelowCutOff(cData3<T>& dataSrc, T cutoff, T value);
    void    setValueAboveCutOff(cData3<T>& dataSrc, T cutoff, T value);

    void    threshold(      cData3<T>& dataSrc,                     T cutoff);
    void    threshold(const cData3<T>& dataSrc, cData3<T>& dataDst, T cutoff);
};

} // namespace gem

#endif

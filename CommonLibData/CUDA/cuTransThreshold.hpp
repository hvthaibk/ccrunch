/***********************************************************************
 *  File:       cuTransThreshold.hpp
 *
 *  Purpose:    Header file for a thresholding transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUTRANS_THRESHOLD_HPP__
#define __GEM_CUTRANS_THRESHOLD_HPP__

#include "cuData3.hpp"

namespace gem {

template <typename T>
class cuTransThreshold
{
private:

public:
     cuTransThreshold(void) {};
    ~cuTransThreshold()     {};

    // thresholding
    void    threshold(      cuData3<T>& dataSrc,                     T cutoff);
    void    threshold(const cuData3<T>& dataSrc, cuData3<T>& dataDst, T cutoff);
};

} // namespace gem

#endif

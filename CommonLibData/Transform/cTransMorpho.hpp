/***********************************************************************
 *  File:       cTransMorpho.hpp
 *
 *  Purpose:    Header file for a transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CTRANS_MORPHO_HPP__
#define __GEM_CTRANS_MORPHO_HPP__

#include "cData3.hpp"

namespace gem {

template <typename T>
class cTransMorpho
{
private:

public:
     cTransMorpho(void) {};
    ~cTransMorpho()     {};

    // dilation
    void    dilate(      cData3<T>& dataSrc,                     T radius);
    void    dilate(const cData3<T>& dataSrc, cData3<T>& dataDst, T radius);
};

} // namespace gem

#endif

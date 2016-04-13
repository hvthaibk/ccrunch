/***********************************************************************
 *  File:       cTIFF.hpp
 *
 *  Purpose:    Header file for a TIFF class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CTIFF_HPP__
#define __GEM_CTIFF_HPP__

#include "cArray3.hpp"

#define TIFF_SUCCESS    1

namespace gem {

template <typename T>
class cTIFF
{
protected:

public:
     cTIFF(void) {};
    ~cTIFF()     {};

    // check extension
    void checkExtension(const std::string& fileName) const;

    // info
    void getInfo(const std::string& fileName) const;

    // input & output
    void read (const std::string& fileName, cArray3<T>& object) const;
    void write(const cArray3<T>& object, const std::string& fileName, eNorm8bit norm8bit = NORM8BIT_FALSE) const;
};

} // namespace gem

#endif

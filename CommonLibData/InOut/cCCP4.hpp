/***********************************************************************
 *  File:       cCCP4.hpp
 *
 *  Purpose:    Header file for a CCP4 class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CCCP4_HPP__
#define __GEM_CCCP4_HPP__

#include "cData3.hpp"
#include "cMapHeader.hpp"

namespace gem {

template <typename T>
class cCCP4
{
protected:

public:
     cCCP4(void) {};
    ~cCCP4()     {};

    // check extension
    void checkExtension(const std::string& fileName) const;

    // header
    void readHeader(const std::string& fileName, cMapHeader& header);

    // data
    void readData(const std::string& fileName, cData3<T>& object, cMapHeader& header);

    // input & output
    void read (const std::string& fileName, cData3<T>& object, cMapHeader* header = NULL);
    void write(const cData3<T>& object, const std::string& fileName, cMapHeader* header = NULL);
};

} // namespace gem

#endif

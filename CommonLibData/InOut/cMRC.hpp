/***********************************************************************
 *  File:       cMRC.hpp
 *
 *  Purpose:    Header file for an MRC class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CMRC_HPP__
#define __GEM_CMRC_HPP__

#include "cData3.hpp"
#include "cMapHeader.hpp"

namespace gem {

template <typename T>
class cMRC
{
protected:

public:
     cMRC(void) {};
    ~cMRC()     {};

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

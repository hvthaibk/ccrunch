/***********************************************************************
 *  File:       cResultFitter.hpp
 *
 *  Purpose:    Header file for gEMfitter's result structure
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CRESULT_FITTER_HPP__
#define __GEM_CRESULT_FITTER_HPP__

#include "cVector3.hpp"

#include <string>
#include <vector>

namespace gem {

template <typename T>
struct sFitterResult
{
    std::string     _numStr;
    cVector3<T>     _pos;
    cVector3<T>     _angle;
    T               _corrVal;

    sFitterResult(void) {};
    sFitterResult(std::string numStr, cVector3<T> pos, cVector3<T> angle, T corrVal) {
        _numStr  = numStr;
        _pos     = pos;
        _angle   = angle;
        _corrVal = corrVal;
    }

    template <typename U>
    friend  std::ostream& operator<<(std::ostream& os, const sFitterResult<U>& result);
};

template <typename T>
class cResultFitter : public std::vector< sFitterResult<T> >
{
protected:

public:
    // operator overloading
    template <typename U>
    friend  std::ostream& operator<<(std::ostream& os, const cResultFitter<U>& result);

    // input & output
    void    read (const std::string& fileName);
    void    write(const std::string& fileName) const;
};

} // namespace gem

#endif

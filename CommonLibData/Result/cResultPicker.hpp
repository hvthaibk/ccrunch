/***********************************************************************
 *  File:       cResultPicker.hpp
 *
 *  Purpose:    Header file for gEMpicker's result structure
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CRESULT_PICKER_HPP__
#define __GEM_CRESULT_PICKER_HPP__

#include <string>
#include <vector>

namespace gem {

template <typename T>
struct sPickerResult
{
    std::string     _numStr;
    size_t          _x, _y;
    size_t          _boxSize;
    size_t          _tplIdx;
    std::string     _tplName, _mskName;
    T               _angle2D;
    T               _corrVal;

    sPickerResult(void) {};
    sPickerResult(std::string numStr, size_t x, size_t y,
                  size_t boxSize, size_t tplIdx,
                  std::string tplName, std::string mskName,
                  T angle2D, T corrVal) {
        _numStr  = numStr;
        _x       = x;
        _y       = y;
        _boxSize = boxSize;
        _tplIdx  = tplIdx;
        _tplName = tplName;
        _mskName = mskName;
        _angle2D = angle2D;
        _corrVal = corrVal;
    }

    template <typename U>
    friend  std::ostream& operator<<(std::ostream& os, const sPickerResult<U>& result);
};

template <typename T>
class cResultPicker : public std::vector< sPickerResult<T> >
{
protected:

public:    
    // operator overloading
    template <typename U>
    friend  std::ostream& operator<<(std::ostream& os, const cResultPicker<U>& result);

    // input & output
    void    read        (const std::string& fileName);
    void    write       (const std::string& fileName) const;
    void    writeBoxEMAN(const std::string& fileName, size_t refNrow) const;

    // remove close particles
    size_t  refineUsingDistance(size_t boxDist);
    // remove particles with too high CC value
    size_t  refineUsingHighThreshold(float thresholdHigh);
};

} // namespace gem

#endif

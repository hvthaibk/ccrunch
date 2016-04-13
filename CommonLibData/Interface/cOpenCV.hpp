/***********************************************************************
 *  File:       cOpenCV.hpp
 *
 *  Purpose:    Header file for an OpenCV class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_COPENCV_HPP__
#define __GEM_COPENCV_HPP__

#include "cData3.hpp"

#ifdef __GEM_USE_OPENCV__

namespace gem {

template <typename T>
class cOpenCV
{
protected:

public:
     cOpenCV(void) {};
    ~cOpenCV()     {};

    // put text onto image
    void putText(const cData3<T>& object, const std::string& fileName,
                 const std::string& textStr,
                 eColorName textColor = COLOR_BLUE,
                 T scale = 1,
                 eNorm8bit norm8bit = NORM8BIT_FALSE);

    // input & output
    void read (const std::string& fileName, cData3<T>& object) const;
    void write(const cData3<T>& object, const std::string& fileName, eNorm8bit norm8bit = NORM8BIT_FALSE) const;

    // viewing
    void view2D(const cData3<T>& object, const std::string& title = "OpenCV view2D", eNorm8bit norm8bit = NORM8BIT_TRUE);
    void view3D(const cData3<T>& object, const std::string& title = "OpenCV view3D");

    void view2D(const std::string& title, unsigned int nImg, unsigned int sep, const cData3<T>* object, ...);
};

} // namespace gem

#endif  // __GEM_USE_OPENCV__

#endif

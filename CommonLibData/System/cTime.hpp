/***********************************************************************
 *  File:       cTime.hpp
 *
 *  Purpose:    Header file for a high-resolution timer
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CTIME_HPP__
#define __GEM_CTIME_HPP__

#include <iostream>
#include <string>
#include <sys/time.h>

namespace gem {

class cTime
{
private:
    timeval       time;
    double        timeStart, timeStop;

public:
     cTime(void);
    ~cTime();

    void tic();
    void toc(std::string str = "");

    double getElapsedTime();
};

} // namespace gem

#endif

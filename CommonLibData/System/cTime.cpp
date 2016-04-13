/***********************************************************************
 *  File:       cTime.cpp
 *
 *  Purpose:    Implementation of a high-resolution timer
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cTime.hpp"

namespace gem {

cTime::cTime(void)
{
    timeStart = 0;
    timeStop  = 0;
}

cTime::~cTime()
{
}

void cTime::tic()
{
    gettimeofday(&time, NULL);

    // convert from microseconds to seconds
    timeStart = (double) time.tv_sec + (double) time.tv_usec/1000000.0;
}

void cTime::toc(std::string str)
{
    gettimeofday(&time, NULL);

    // convert from microseconds to seconds
    timeStop = (double) time.tv_sec + (double) time.tv_usec/1000000.0;

    if (str.empty()) {
        std::cout << "     elapsed time: "
                  << timeStop - timeStart << std::endl;
    }
    else if (str.compare("noshow") == 0) {
    }
    else {
        std::cout << str << timeStop - timeStart << std::endl;
    }
}

double cTime::getElapsedTime()
{
    return timeStop - timeStart;
}

} // namespace gem

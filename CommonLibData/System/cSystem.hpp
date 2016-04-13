/***********************************************************************
 *  File:       cSystem.hpp
 *
 *  Purpose:    Header file for a system-info class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CSYSTEM_HPP__
#define __GEM_CSYSTEM_HPP__

#include "config.hpp"
#include "macro.hpp"

#include <iostream>
#include <cstdio>
#include <cstring>

namespace gem {

class cSystem
{
private:

public:
     cSystem(void) {};
    ~cSystem()     {};

    const std::string   getComputerName  (void);
    unsigned int        getNumProcessors (void);
    unsigned int        getNumCoresPerCPU(void);
    unsigned int        getNumGPUs       (void);
    bool                getEndianness    (void);
    void                getCacheInfo     (void);

    void    startOpenMP(void);
    void    stopOpenMP (void);

    void    startNthreadsFFTW(void);
    void    stopNthreadsFFTW (void);

    void    checkComputingResource(unsigned int nGPU, unsigned int nCPU);
    void    checkMemoryLeak       (void);
};

} // namespace gem

#endif

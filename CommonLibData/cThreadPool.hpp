/***********************************************************************
 *  File:       cThreadPool.hpp
 *
 *  Purpose:    Header file for thread-safe data structures
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CTHREADPOOL_HPP__
#define __GEM_CTHREADPOOL_HPP__

#include "config.hpp"
#include "macro.hpp"

#include "boost/thread.hpp"
#include <queue>

namespace gem {

template<class T>
class cThreadPool
{
private:

public:
     cThreadPool(void) {}
    ~cThreadPool()     {}
};

} // namespace gem

#endif

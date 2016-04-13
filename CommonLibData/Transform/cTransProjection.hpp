/***********************************************************************
 *  File:       cTransProjection.hpp
 *
 *  Purpose:    Header file for a 2D/3D projection class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CTRANS_PROJECTION_HPP__
#define __GEM_CTRANS_PROJECTION_HPP__

#include "cData3.hpp"

namespace gem {

template <typename T>
class cTransProjection
{
private:

public:
     cTransProjection(void) {};
    ~cTransProjection()     {};

    // projection
    void    project2D(const cData3<T>& dataSrc, cData3<T>& dataDst,
                      std::vector<T>& theta,
                      eInter inter = INTER_LINEAR);
    void    project3D(const cData3<T>& dataSrc, cData3<T>& dataDst,
                      std::vector<T>& alpha,
                      std::vector<T>& beta,
                      eInter inter = INTER_LINEAR);

    // tracing
    void    trace2D(const cData3<T>& dataSrc, cData3<T>& dataDst,
                    std::vector<T>& theta,
                    eTrace trace = TRACE_SUM,
                    eInter inter = INTER_LINEAR);
    void    trace3D(const cData3<T>& dataSrc, cData3<T>& dataDst,
                    std::vector<T>& alpha,
                    std::vector<T>& beta,
                    std::vector<T>& gamma,
                    eTrace trace = TRACE_SUM,
                    eInter inter = INTER_LINEAR);
};

} // namespace gem

#endif

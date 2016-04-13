/***********************************************************************
 *  File:       cuTransSimilarity.hpp
 *
 *  Purpose:    Header file for a similarity transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUTRANS_SIMILARITY_HPP__
#define __GEM_CUTRANS_SIMILARITY_HPP__

#include "cuData3.hpp"

namespace gem {

template <typename T>
class cuTransSimilarity
{
private:

public:
     cuTransSimilarity(void) {};
    ~cuTransSimilarity()     {};

    void    extendRotate(      cuData3<T>& dataSrc,                      eExpRot exprot = EXP_ROT_VALID);
    void    extendRotate(const cuData3<T>& dataSrc, cuData3<T>& dataDst, eExpRot exprot = EXP_ROT_VALID);

    void    translate   (      cuData3<T>& dataSrc,                      const cVector3<T>& offset, eInter inter = INTER_LINEAR);
    void    translate   (const cuData3<T>& dataSrc, cuData3<T>& dataDst, const cVector3<T>& offset, eInter inter = INTER_LINEAR);

    void    rotate      (      cuData3<T>& dataSrc,                      const cVector3<T>& angle,  eInter inter = INTER_LINEAR);
    void    rotate      (const cuData3<T>& dataSrc, cuData3<T>& dataDst, const cVector3<T>& angle,  eInter inter = INTER_LINEAR);

    void    scale       (      cuData3<T>& dataSrc,                      T factor,                  eInter inter = INTER_LINEAR);
    void    scale       (const cuData3<T>& dataSrc, cuData3<T>& dataDst, T factor,                  eInter inter = INTER_LINEAR);

    void    rotTrans    (      cuData3<T>& dataSrc,                      const cVector3<T>& angle, const cVector3<T>& offset, eInter inter = INTER_LINEAR);
    void    rotTrans    (const cuData3<T>& dataSrc, cuData3<T>& dataDst, const cVector3<T>& angle, const cVector3<T>& offset, eInter inter = INTER_LINEAR);

    void    interp(const cuData3<T>& R,
                   const cuData3<T>& dataSrc, cuData3<T>& dataDst,
                   eInter inter = INTER_LINEAR);
    void    interp(const cuData3<T>& R, const cuData3<T>& C,
                   const cuData3<T>& dataSrc, cuData3<T>& dataDst,
                   eInter inter = INTER_LINEAR);
    void    interp(const cuData3<T>& R, const cuData3<T>& C, const cuData3<T>& S,
                   const cuData3<T>& dataSrc, cuData3<T>& dataDst,
                   eInter inter = INTER_LINEAR);
};

} // namespace gem

#endif

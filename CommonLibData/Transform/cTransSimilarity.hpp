/***********************************************************************
 *  File:       cTransSimilarity.hpp
 *
 *  Purpose:    Header file for a similarity transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CTRANS_SIMILARITY_HPP__
#define __GEM_CTRANS_SIMILARITY_HPP__

#include "cData3.hpp"

namespace gem {

template <typename T>
class cTransSimilarity
{
private:

public:
     cTransSimilarity(void) {};
    ~cTransSimilarity()     {};

    void    bin         (      cData3<T>& dataSrc,                     size_t factor);
    void    bin         (const cData3<T>& dataSrc, cData3<T>& dataDst, size_t factor);

    void    extendRotate(      cData3<T>& dataSrc,                     eExpRot exprot = EXP_ROT_VALID);
    void    extendRotate(const cData3<T>& dataSrc, cData3<T>& dataDst, eExpRot exprot = EXP_ROT_VALID);

    void    translate   (      cData3<T>& dataSrc,                     const cVector3<T>& offset, eInter inter = INTER_LINEAR);
    void    translate   (const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& offset, eInter inter = INTER_LINEAR);

    void    rotate      (      cData3<T>& dataSrc,                     const cVector3<T>& angle,  eInter inter = INTER_LINEAR);
    void    rotate      (const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& angle,  eInter inter = INTER_LINEAR);

    void    scale       (      cData3<T>& dataSrc,                     T factor, eInter inter = INTER_LINEAR);
    void    scale       (const cData3<T>& dataSrc, cData3<T>& dataDst, T factor, eInter inter = INTER_LINEAR);
    void    scale       (      cData3<T>& dataSrc,                     const cVector3<T>& origin, T factor, eInter inter = INTER_LINEAR);
    void    scale       (const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& origin, T factor, eInter inter = INTER_LINEAR);

    void    rotTrans    (      cData3<T>& dataSrc,                     const cVector3<T>& angle, const cVector3<T>& offset, eInter inter = INTER_LINEAR);
    void    rotTrans    (const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& angle, const cVector3<T>& offset, eInter inter = INTER_LINEAR);

    void    interp(const cData3<T>& R,
                   const cData3<T>& dataSrc, cData3<T>& dataDst,
                   eInter inter = INTER_LINEAR);
    void    interp(const cData3<T>& R, const cData3<T>& C,
                   const cData3<T>& dataSrc, cData3<T>& dataDst,
                   eInter inter = INTER_LINEAR);
    void    interp(const cData3<T>& R, const cData3<T>& C, const cData3<T>& S,
                   const cData3<T>& dataSrc, cData3<T>& dataDst,
                   eInter inter = INTER_LINEAR);
};

} // namespace gem

#endif

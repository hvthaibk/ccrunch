/***********************************************************************
 *  File:       cStackRegTrans.hpp
 *
 *  Purpose:    Header file for a stack registration class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CSTACK_REG_TRANS_HPP__
#define __GEM_CSTACK_REG_TRANS_HPP__

#include "xcorr.hpp"
#include "cData3.hpp"

namespace gem {

template <typename T>
class cStackRegTrans
{
private:
    std::vector< cVector3<T> >  _corrRes;

    cSize3                      _tplSize, _tplOffset;
    cSize3                      _refSize, _refOffset;

public:
     cStackRegTrans(void) {};
    ~cStackRegTrans()     {};

    void    memFree(void);

    std::vector< cVector3<T> >&  getResult(void);
    void                         setResult(const std::vector< cVector3<T> >& corrRes);

    void        setROIs(const cSize3& tplSize, const cSize3& tplOffset,
                        const cSize3& refSize, const cSize3& refOffset);

    // translational registration
    void    computeOffsetsUsingCorr(const std::string& dirDataIn,
                                    const std::vector<std::string>& fileList,
                                    eCCmode mode = CC_MODE_NCORR);
    void    alignSections(const std::string& dirDataIn,
                          const std::string& dirDataOut,
                          const std::vector<std::string>& fileList,
                          bool  extend = true);
};

} // namespace gem

#endif

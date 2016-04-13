/***********************************************************************
 *  File:       cStackRegTPS2.hpp
 *
 *  Purpose:    Header file for a stack registration class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CSTACK_REG_TPS2_HPP__
#define __GEM_CSTACK_REG_TPS2_HPP__

#include "xcorr.hpp"
#include "cData3.hpp"

namespace gem {

template <typename T>
class cStackRegTPS2
{
private:
    cSize3      _tplSizeHalf, _searchRange, _step, _border;

public:
     cStackRegTPS2(void) {};
    ~cStackRegTPS2()     {};

    // params
    void    setParamsPointTracking(const cSize3& tplSizeHalf, const cSize3& searchRange,
                                   const cSize3& step,        const cSize3& border);
    void    setParams(void);

    // operations
    void    trackPointsUsingCorr(const std::string& dirDataIn,
                                 const std::string& dirDataTmp,
                                 const std::vector<std::string>& fileList,
                                 eCCmode mode = CC_MODE_NCORR);

    void    computeMatrixW(const std::string& dirDataTmp,
                           const std::vector<std::string>& fileList);

    void    computeMappingAndInterp(const std::string& dirDataIn,
                                    const std::string& dirDataTmp,
                                    const std::string& dirDataOut,
                                    const std::vector<std::string>& fileList);
};

} // namespace gem

#endif

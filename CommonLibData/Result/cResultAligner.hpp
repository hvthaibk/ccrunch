/***********************************************************************
 *  File:       cResultAligner.hpp
 *
 *  Purpose:    Header file for gEMaligner's result structure
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CRESULT_ALIGNER_HPP__
#define __GEM_CRESULT_ALIGNER_HPP__

#include "cData3.hpp"
#include "cSetPoint3.hpp"

#include <string>

namespace gem {

template <typename T>
class cResultAligner
{
protected:

public:
    // convert from row-col to X-Y for TPS2 computation
    void    swapPointsCoordinates(cData3<T>& points);

    // input & output
    void     readCorrVals(const std::string& fileName, cData3<T>& corrs) const;
    void    writeCorrVals(const cData3<T>& corrs, const std::string& fileName) const;

    void     readPoints(const std::string& fileName, cData3<T>& points) const;
    void    writePoints(const cData3<T>& points, const std::string& fileName, unsigned int precision = 4) const;

    // rigid alignment
    void     readResultRigid(const std::string& fileName, std::vector< cVector3<T> >& result) const;
    void    writeResultRigid(const std::vector< cVector3<T> >& result, const std::string& fileName, unsigned int precision = 4) const;
    void    checkResultRigid(const std::vector< cVector3<T> >& result, const cSize3& imSize, T thrOffset, T thrCorrVal) const;
};

} // namespace gem

#endif

/***********************************************************************
 *  File:       cESBTL.hpp
 *
 *  Purpose:    Header file for an ESBTL class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CESBTL_HPP__
#define __GEM_CESBTL_HPP__

#include "cSetPoint3.hpp"

#include <string>

namespace gem {

class cESBTL
{
protected:

public:
     cESBTL(void) {};
    ~cESBTL()     {};

    // check extension
    void   checkExtension(const std::string& fileName) const;

    // input & output
    void   pdbPrintAtom   (const std::string& fileIn, const std::string& atomName = "");
    void   pdbReadAtom    (const std::string& fileIn, cSetPoint3<double>& atomList, const std::string& atomName = "CA");
    void   pdbExtractChain(const std::string& fileIn, const std::string& fileOut, const std::string& chainIDs);

    // transform
    void   pdbRotate       (const std::string& fileIn, const std::string& fileOut, const cVector3<double>& angle);
    void   pdbTranslate    (const std::string& fileIn, const std::string& fileOut, const cVector3<double>& offset);
    void   pdbTransRotTrans(const std::string& fileIn, const std::string& fileOut, const cVector3<double>& offset1,
                                                                                   const cVector3<double>& angle,
                                                                                   const cVector3<double>& offset2);

    // statistics
    double pdbComputeRMSD  (const std::string& fileIn1, const std::string& fileIn2, const std::string& atomName = "CA");
    void   pdbComputeSD    (const std::string& fileIn1, const std::string& fileIn2, double& sqrDist, size_t& numAtom, const std::string& atomName = "CA");
};

} // namespace gem

#endif

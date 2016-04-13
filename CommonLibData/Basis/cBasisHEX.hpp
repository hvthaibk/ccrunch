/***********************************************************************
 *  File:       cBasisHEX.hpp
 *
 *  Purpose:    Header file for a basis-related class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CBASIS_HEX_HPP__
#define __GEM_CBASIS_HEX_HPP__

#include "cData3.hpp"

#ifdef __GEM_USE_HEX__
#include "HEX/hex_molecule.h"
#endif

namespace gem {

template <typename T>
class cBasisHEX
{
private:

public:
     cBasisHEX(void) {};
    ~cBasisHEX()     {};

#ifdef __GEM_USE_HEX__
    void    hexComputeCoeff(const cVector3<double>& origin,
                            const cVector3<double>& spacing,
                            unsigned int nmax, double scale,
                            cData3<double>& coff) const;
    void    hexReconData   (const cVector3<double>& origin,
                            const cVector3<double>& spacing,
                            unsigned int nmax, double scale,
                            const cData3<double>& coff);

    void    hexSetupOption(DockSpec& ds);
    void    hexMatcher    (DockSpec& ds, const cBasisHEX<T>& other, Docking& matching);
#endif
};

} // namespace gem

#endif

/***********************************************************************
 *  File:       cBasisHermite.hpp
 *
 *  Purpose:    Header file for a basis-related class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CBASIS_HERMITE_HPP__
#define __GEM_CBASIS_HERMITE_HPP__

#include "cBasis.hpp"

namespace gem {

template <typename T>
class cBasisHermite : public cBasis<T>
{
private:
    cData3<T>       _basisRow, _basisCol, _basisSec;

public:
     cBasisHermite(void);
    ~cBasisHermite()    ;

    void    memFree(void);

    const cData3<T>& getBasisRow(void) { return _basisRow; }
    const cData3<T>& getBasisCol(void) { return _basisCol; }
    const cData3<T>& getBasisSec(void) { return _basisSec; }

    // debugging
    void    printData(std::string message = "", int width = 8, int precision = 4, bool scientific = 0) const;
    void    printSize(std::string message = "", int width = 8, int precision = 4) const;

    // generate basis functions
    void    genBasisRow(T sigma, eBasisAccuracy accuracy = BASIS_ACCU_MP);
    void    genBasisCol(T sigma, eBasisAccuracy accuracy = BASIS_ACCU_MP);
    void    genBasisSec(T sigma, eBasisAccuracy accuracy = BASIS_ACCU_MP);

    void    genBasis2D(const cInt3& order, cData3<T>& basis2D);
    void    genBasis3D(const cInt3& order, cData3<T>& basis3D);

    // decomposition / reconstruction
    void    computeCoef(const cData3<T>&          data,
                        const std::vector<cInt3>& order,
                              cData3<T>&          coefs);
    void    reconData  (const cData3<T>&          coefs,
                        const std::vector<cInt3>& order,
                              cData3<T>&          data);
};

} // namespace gem

#endif

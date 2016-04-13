/***********************************************************************
 *  File:       cArmadillo.hpp
 *
 *  Purpose:    Header file for an Armadillo-interface class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CARMADILLO_HPP__
#define __GEM_CARMADILLO_HPP__

#include "cData3.hpp"

#ifdef __GEM_USE_ARMADILLO__

#define ARMA_64BIT_WORD
#define ARMA_USE_HDF5

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wcast-qual"
#include "armadillo"
#pragma GCC diagnostic pop

namespace gem {

template <typename T>
class cArmadillo
{
private:

public:
     cArmadillo(void) {};
    ~cArmadillo()     {};

    // check extension
    void checkExtension(const std::string& fileName) const;

    // data conversion
    void    convertData(const cData3<T>& a, arma::Mat<T>& b);
    void    convertData(const arma::Mat<T>& a, cData3<T>& b);

    // input / output
    void    load(const std::string& fileName, cData3<T>& data);
    void    save(const cData3<T>& data, const std::string& fileName);

    // linear algebra
    void    inverseMat (const cData3<T>& A, cData3<T>& B);
    void    multiplyMat(const cData3<T>& A, const cData3<T>& B, cData3<T>& C);
    void    solveAxb   (const cData3<T>& A, const cData3<T>& B, cData3<T>& X);
};

} // namespace gem

#endif  // __GEM_USE_ARMADILLO__

#endif

/***********************************************************************
 *  File:       cEigen.hpp
 *
 *  Purpose:    Header file for an Eigen-interface class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CEIGEN_HPP__
#define __GEM_CEIGEN_HPP__

#include "cData3.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <Eigen/Dense>
#pragma GCC diagnostic pop

namespace gem {

enum eEigenSolverAxb
{
    EIGEN_SOLVER_AXB_LDLT                 = 0,
    EIGEN_SOLVER_AXB_LLT                  = 1,
    EIGEN_SOLVER_AXB_LU                   = 2,
    EIGEN_SOLVER_AXB_FULLPIVLU            = 3,
    EIGEN_SOLVER_AXB_COLPIVHOUSEHOLDERQR  = 4,
    EIGEN_SOLVER_AXB_FULLPIVHOUSEHOLDERQR = 5,
    EIGEN_SOLVER_AXB_HOUSEHOLDERQR        = 6
};

template <typename T>
class cEigen
{
private:

public:
     cEigen(void) {};
    ~cEigen()     {};

    // data conversion
    void    convertData(const cData3<T>& a, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& b);
    void    convertData(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& a, cData3<T>& b);

    void    solveAxb(const cData3<T>& A, const cData3<T>& B, cData3<T>& X,
                     eEigenSolverAxb solver = EIGEN_SOLVER_AXB_LU);
};

} // namespace gem

#endif

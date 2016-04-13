/***********************************************************************
 *  File:       cEigen.cpp
 *
 *  Purpose:    Implementation of an Eigen-interface class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cEigen.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
#pragma GCC diagnostic pop

namespace gem {

template <typename T>
void cEigen<T>::convertData(const cData3<T>& a,
                            Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& b)
{
    const std::string   funcName("void cEigen<T>::convertData("
                                    "const cData3<T>& a, "
                                    "Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& b)");

    a.requireNonEmpty(funcName);
    require(a.getNsec() == 1, funcName + ": a is not a matrix");

    b.resize(a.getNrow(),a.getNcol());

    for (size_t i = 0, k = 0; i < a.getNrow(); i++) {
        for (size_t j = 0; j < a.getNcol(); j++, k++) {
            b(i,j) = a[k];
        }
    }
}

template <typename T>
void cEigen<T>::convertData(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& a,
                            cData3<T>& b)
{
    const std::string   funcName("void cEigen<T>::convertData("
                                    "const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& a, "
                                    "cData3<T>& b)");

    require(a.size() > 0, funcName + ": a is not a matrix");

    b.memReAlloc(cSize3(a.rows(),a.cols(),1));

    for (size_t i = 0, k = 0; i < b.getNrow(); i++) {
        for (size_t j = 0; j < b.getNcol(); j++, k++) {
            b[k] = a(i,j);
        }
    }
}

template <typename T>
void cEigen<T>::solveAxb(const cData3<T>& A, const cData3<T>& B, cData3<T>& X,
                         eEigenSolverAxb solver)
{
    const std::string   funcName("void cEigen::solveAxb("
                                    "const cData3<T>& A, "
                                    "const cData3<T>& B, "
                                    "cData3<T>& X, "
                                    "eEigenSolverAxb solver)");

    require(A.getDimension() == 2 && A.getNsec() == 1,
            funcName + ": A is not a matrix");
    require(A.getNcol() == B.getNrow(),
            funcName + ": A and B are not compatible");

    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>  _A, _B, _X;

    convertData(A,_A);
    convertData(B,_B);

    switch (solver) {
        case EIGEN_SOLVER_AXB_LDLT:
            _X = _A.ldlt().solve(_B);
            break;
        case EIGEN_SOLVER_AXB_LLT:
            _X = _A.llt().solve(_B);
            break;
        case EIGEN_SOLVER_AXB_LU:
            _X = _A.lu().solve(_B);
            break;
        case EIGEN_SOLVER_AXB_FULLPIVLU:
            _X = _A.fullPivLu().solve(_B);
            break;
        case EIGEN_SOLVER_AXB_COLPIVHOUSEHOLDERQR:
            _X = _A.colPivHouseholderQr().solve(_B);
            break;
        case EIGEN_SOLVER_AXB_FULLPIVHOUSEHOLDERQR:
            _X = _A.fullPivHouseholderQr().solve(_B);
            break;
        case EIGEN_SOLVER_AXB_HOUSEHOLDERQR:
            _X = _A.householderQr().solve(_B);
            break;
        default:
            ERROR(funcName, "unsupported solver");
    }

    convertData(_X,X);
}

// instantiation
template class cEigen<float >;
template class cEigen<double>;

} // namespace gem

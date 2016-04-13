/***********************************************************************
 *  File:       cPCA.hpp
 *
 *  Purpose:    Header file for PCA classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CPCA_HPP__
#define __GEM_CPCA_HPP__

#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#include <Eigen/Dense>
#pragma GCC diagnostic pop

namespace gem {

class cPCA2
{
private:
    std::vector<float>  _x;         // Initial matrix as vector filled by row
    Eigen::MatrixXf     _xXf;       // Initial matrix as Eigen MatrixXf structure

    size_t              _nrow,      // Number of rows in matrix x
                        _ncol;      // Number of cols in matrix x

    bool                _is_center, // Whether the variables should be shifted to be zero centered
                        _is_scale,  // Whether the variables should be scaled to have unit variance
                        _is_corr;   // PCA with correlation matrix, not covariance

    std::string         _method;    // svd, cor, cov
                        // There are different methods used. The most used is SVD.
                        // But in some cases it may be correlation or covariance matrices.
                        //      "svd" - PCA with singular value decomposition
                        //      "cor" - PCA with correlation matrix
                        //      "cov" - PCA with covariance matrix

    std::vector<size_t> _eliminated_columns;    // Numbers of eliminated columns
                        // If standard deviation of a column is equal to 0,
                        // the column shoud be rejected, or PCA will fail

    std::vector<float>  _sd,            // Standard deviation of each component
                        _prop_of_var,   // Proportion of variance of each component
                        _cum_prop;      // Cumulative proportion

    std::vector<float>  _scores;        // Rotated values
                        // Calculated scores (coordinates in a new space)
                        // as vector. Matrix filled by row

    unsigned int        _kaiser,        // Number of PC according Kaiser criterion
                        _thresh95;      // Number of PC according 95% variance threshold

public:
     cPCA2(void);
    ~cPCA2()    ;

    size_t  getNrow(void) const { return _nrow; }
    size_t  getNcol(void) const { return _ncol; }

    bool    isCenter(void) const { return _is_center; }
    bool    isScale (void) const { return _is_scale;  }
    bool    isCorr  (void) const { return _is_corr;   }

    std::string         getMethodName    (void) const { return _method; }
    std::vector<size_t> getEliminatedCols(void) const { return _eliminated_columns; }

    std::vector<float>  getSd     (void) const { return _sd;         };
    std::vector<float>  getVarProp(void) const {return _prop_of_var; };
    std::vector<float>  getVarCum (void) const { return _cum_prop;   };
    std::vector<float>  getScores (void) const { return _scores;     };

    unsigned int        getNPCkaiser  (void) const { return _kaiser;   };
    unsigned int        getNPCthresh95(void) const { return _thresh95; };

    /*!
        The main method for performin Principal Component Analysis
        \param  x           Initial data matrix
        \param  nrow        Number of matrix row
        \param  ncol        Number of matrix col
        \param  is_center   Whether the variables should be shifted to be zero centered
        \param  is_scale    Whether the variables should be scaled to have unit variance
        \param  is_corr     Correlation matrix will be used instead of covariance matrix
        \result
            0  if everything is Ok
            -1 if there were some errors
    */
    int computePC(std::vector<float>& x,
                  size_t nrow,
                  size_t ncol,
                  bool is_center = true,
                  bool is_scale  = true,
                  bool is_corr   = true);
};

class cPCA3
{
private:
    std::vector<float>  _x;         // Initial matrix as vector filled by row
    Eigen::MatrixXf     _xXf;       // Initial matrix as Eigen MatrixXf structure

    size_t              _nrow,      // Number of row in matrix x
                        _ncol,      // Number of cols in matrix x
                        _nsec;      // Number of secs in matrix x

    bool                _is_center, // Whether the variables should be shifted to be zero centered
                        _is_scale,  // Whether the variables should be scaled to have unit variance
                        _is_corr;   // PCA with correlation matrix, not covariance

    std::string         _method;    // svd, cor, cov
                        // There are different methods used. The most used is SVD.
                        // But in some cases it may be correlation or covariance matrices.
                        //      "svd" - PCA with singular value decomposition
                        //      "cor" - PCA with correlation matrix
                        //      "cov" - PCA with covariance matrix

    std::vector<size_t> _eliminated_columns;    // Numbers of eliminated columns
                        // If standard deviation of a column is equal to 0,
                        // the column shoud be rejected, or PCA will fail

    std::vector<float>  _sd,            // Standard deviation of each component
                        _prop_of_var,   // Proportion of variance of each component
                        _cum_prop;      // Cumulative proportion

    std::vector<float>  _scores;        // Rotated values
                        // Calculated scores (coordinates in a new space)
                        // as vector. Matrix filled by row

    unsigned int        _kaiser,        // Number of PC according Kaiser criterion
                        _thresh95;      // Number of PC according 95% variance threshold

public:
     cPCA3(void);
    ~cPCA3()    ;

    size_t  getNrow(void) const { return _nrow; }
    size_t  getNcol(void) const { return _ncol; }
    size_t  getNsec(void) const { return _nsec; }

    bool    isCenter(void) const { return _is_center; }
    bool    isScale (void) const { return _is_scale;  }
    bool    isCorr  (void) const { return _is_corr;   }

    std::string         getMethodName    (void) const { return _method; }
    std::vector<size_t> getEliminatedCols(void) const { return _eliminated_columns; }

    std::vector<float>  getSd     (void) const { return _sd;         };
    std::vector<float>  getVarProp(void) const {return _prop_of_var; };
    std::vector<float>  getVarCum (void) const { return _cum_prop;   };
    std::vector<float>  getScores (void) const { return _scores;     };

    unsigned int        getNPCkaiser  (void) const { return _kaiser;   };
    unsigned int        getNPCthresh95(void) const { return _thresh95; };

    /*!
        The main method for performin Principal Component Analysis
        \param  x           Initial data matrix
        \param  nrow        Number of matrix row
        \param  ncol        Number of matrix col
        \param  nsec        Number of matrix sec
        \param  is_center   Whether the variables should be shifted to be zero centered
        \param  is_scale    Whether the variables should be scaled to have unit variance
        \param  is_corr     Correlation matrix will be used instead of covariance matrix
        \result
            0  if everything is Ok
            -1 if there were some errors
    */
    int computePC(std::vector<float>& x,
                  size_t nrow,
                  size_t ncol,
                  size_t nsec,
                  bool is_center = true,
                  bool is_scale  = true,
                  bool is_corr   = true);
};

} // namespace gem

#endif

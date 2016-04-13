/***********************************************************************
 *  File:       transform_label.cpp
 *
 *  Purpose:    Implementation of transformation functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "transform.hpp"

namespace gem {

/* REFERENCE
 * http://www.mathworks.com/matlabcentral/fileexchange/26946-label-connected-components-in-2-d-array
 * http://www.mathworks.com/matlabcentral/fileexchange/16938-region-adjacency-graph-rag
 * Union-Find Algorithms: http://www.cs.duke.edu/courses/cps100e/fall09/notes/UnionFind.pdf
 */

/**************************
 * Labeling
 *************************/

// 1D
template <typename T>
void transform_label(const T*      const arraySrc, size_t nRow,
                           size_t* const labelDst)
{
    assert(arraySrc != NULL);
    assert(labelDst != NULL);
    assert(nRow > 0);
}

// instantiation
template
void transform_label<float >(const float*  const arraySrc, size_t nRow,
                                   size_t* const labelDst);
template
void transform_label<double>(const double* const arraySrc, size_t nRow,
                                   size_t* const labelDst);

// 2D
template <typename T>
void transform_label(const T*      const arraySrc, size_t nRow, size_t nCol,
                           size_t* const labelDst)
{
    assert(arraySrc != NULL);
    assert(labelDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
}

// instantiation
template
void transform_label<float >(const float*  const arraySrc, size_t nRow, size_t nCol,
                                   size_t* const labelDst);
template
void transform_label<double>(const double* const arraySrc, size_t nRow, size_t nCol,
                                   size_t* const labelDst);

// 3D
template <typename T>
void transform_label(const T*      const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                           size_t* const labelDst)
{
    assert(arraySrc != NULL);
    assert(labelDst != NULL);
    assert(nRow > 0);
    assert(nCol > 0);
    assert(nSec > 0);
}

// instantiation
template
void transform_label<float >(const float*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                   size_t* const labelDst);
template
void transform_label<double>(const double* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                   size_t* const labelDst);

} // namespace gem

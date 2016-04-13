/***********************************************************************
 *  File:       transform_watershed.cpp
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

#include <map>

namespace gem {

/**************************
 * Watershed
 *************************/

template <typename T>
static size_t transform_watershed_oldseed(size_t mRow, std::vector<size_t>& seedPointRow,
                                          size_t mCol, std::vector<size_t>& seedPointCol,
                                          size_t mSec, std::vector<size_t>& seedPointSec,
                                          size_t nDim, T seedDistMin)
{
    switch (nDim) {
        case 1:
            for (size_t i = 0; i < seedPointRow.size(); i++) {
                if (dist((T)mRow,(T)seedPointRow[i]) <= seedDistMin) {
                    return i+1;
                }
            }
            break;
        case 2:
            for (size_t i = 0; i < seedPointRow.size(); i++) {
                if (dist((T)mRow,(T)mCol,(T)seedPointRow[i],(T)seedPointCol[i]) <= seedDistMin) {
                    return i+1;
                }
            }
            break;
        case 3:
            for (size_t i = 0; i < seedPointRow.size(); i++) {
                if (dist((T)mRow,(T)mCol,(T)mSec,(T)seedPointRow[i],(T)seedPointCol[i],(T)seedPointSec[i]) <= seedDistMin) {
                    return i+1;
                }
            }
            break;
        default:
            ERROR("transform_watershed_newseed", "unsupported dimension");
    }

    seedPointRow.push_back(mRow);
    seedPointCol.push_back(mCol);
    seedPointSec.push_back(mSec);

    return 0;
}

// 1D/2D/3D
template <typename T>
void transform_watershed(const T*      const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                               size_t* const arrayLbl,
                               size_t        nDim,
                               T             threshold,
                               T             seedDistMin,
                               eWaterDir     water,
                               std::vector< std::set<size_t> >& neighborList)
{
    assert(arraySrc != NULL && arrayLbl != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    typedef typename std::multimap< T, size_t, std::less<T> >   elementType;

    T                       value;
    size_t                  ind, indAdj;
    elementType             elementList;
    size_t                  arrayLength;
    std::vector<size_t>     seedPointRow, seedPointCol, seedPointSec;

    switch (nDim) {
        case 1:
            arrayLength = nRow;
            break;
        case 2:
            arrayLength = nRow * nCol;
            break;
        case 3:
            arrayLength = nRow * nCol * nSec;
            break;
        default:
            ERROR("transform_watershed", "unsupported dimension");
    }

    // prepare all elements having intensity above/below threshold
    switch (water) {
        case WATER_DIR_RISE:
            for (ind = 0; ind < arrayLength; ind++) {
                value = arraySrc[ind];

                if (value < threshold) {
                    elementList.insert(std::pair< T, size_t >(value, ind));
                }
            }
            break;
        case WATER_DIR_FALL:
            for (ind = 0; ind < arrayLength; ind++) {
                value = arraySrc[ind];

                if (value > threshold) {
                    elementList.insert(std::pair< T, size_t >(-value, ind));
                }
            }
            break;
    }

    // traverse the element list in intensity-decreasing order
    size_t                      iLabel = 1;
    std::map<size_t,size_t>     labelAdj;
    size_t                      labelNearBy;
    ptrdiff_t                   mRow = 0, mCol = 0, mSec = 0;
    for (typename elementType::iterator iterElement  = elementList.begin();
                                        iterElement != elementList.end();
                                        iterElement++) {
        value = iterElement->first;
        ind   = iterElement->second;

        labelAdj.clear();

        // neighboring elements
        switch (nDim) {
            case 1:     // 2-connected
                mRow = ind;

                for (ptrdiff_t iRow = mRow-1; iRow <= mRow+1; iRow++) {
                    if (iRow < 0 || iRow >= (ptrdiff_t) nRow) continue;

                    indAdj = iRow;

                    if (arrayLbl[indAdj] > 0) {
                        labelAdj[arrayLbl[indAdj]] += 1;
                    }
                }
                break;
            case 2:     // 8-connected
                mRow = ind / nCol;
                mCol = ind % nCol;

                for (ptrdiff_t iRow = mRow-1; iRow <= mRow+1; iRow++) {
                    if (iRow < 0 || iRow >= (ptrdiff_t) nRow) continue;

                    for (ptrdiff_t iCol = mCol-1; iCol <= mCol+1; iCol++) {
                        if (iCol < 0 || iCol >= (ptrdiff_t) nCol) continue;

                        indAdj = iRow*nCol + iCol;

                        if (arrayLbl[indAdj] > 0) {
                            labelAdj[arrayLbl[indAdj]] += 1;
                        }
                    }
                }
                break;
            case 3:     // 26-connected
                mRow = ind / (nCol*nSec);
                mCol = ind / nSec - mRow*nCol;
                mSec = ind % nSec;

                for (ptrdiff_t iRow = mRow-1; iRow <= mRow+1; iRow++) {
                    if (iRow < 0 || iRow >= (ptrdiff_t) nRow) continue;

                    for (ptrdiff_t iCol = mCol-1; iCol <= mCol+1; iCol++) {
                        if (iCol < 0 || iCol >= (ptrdiff_t) nCol) continue;

                        for (ptrdiff_t iSec = mSec-1; iSec <= mSec+1; iSec++) {
                            if (iSec < 0 || iSec >= (ptrdiff_t) nSec) continue;

                            indAdj = (iRow*nCol + iCol) * nSec + iSec;

                            if (arrayLbl[indAdj] > 0) {
                                labelAdj[arrayLbl[indAdj]] += 1;
                            }
                        }
                    }
                }
                break;
            default:
                ERROR("transform_watershed", "unsupported dimension");
        }

        if (labelAdj.size() == 0) {
        // no adjacent regions found
            labelNearBy = transform_watershed_oldseed(mRow, seedPointRow,
                                                      mCol, seedPointCol,
                                                      mSec, seedPointSec,
                                                      nDim, seedDistMin);

            if (labelNearBy) {
                // the element is added to a nearby region
                arrayLbl[ind] = labelNearBy;
            }
            else {
                // the element starts a new region
                arrayLbl[ind] = iLabel;
                neighborList.resize(iLabel);
                iLabel++;
            }
        }
        else if (labelAdj.size() == 1) {
        // one adjacent region found, the element is added to the region
            arrayLbl[ind] = labelAdj.begin()->first;
        }
        else {
        // two or more adjacent regions found, the element is assigned to the region having the most adjacent elements
            size_t      numAdj = 0;
            size_t      numAdjMaxLabel = 0;

            for (std::map<size_t,size_t>::iterator iterLabel  = labelAdj.begin();
                                                   iterLabel != labelAdj.end();
                                                   iterLabel++) {
                if (iterLabel->second > numAdj) {
                    numAdj         = iterLabel->second;
                    numAdjMaxLabel = iterLabel->first;
                }
            }

            arrayLbl[ind] = numAdjMaxLabel;

            // neighborList
            for (std::map<size_t,size_t>::iterator iterLabel  = labelAdj.begin();
                                                   iterLabel != labelAdj.end();
                                                   iterLabel++) {
                if (iterLabel->first != numAdjMaxLabel) {
                    neighborList[numAdjMaxLabel-1].insert(iterLabel->first);
                    neighborList[iterLabel->first-1].insert(numAdjMaxLabel);
                }
            }
        }
    }
}

// instantiation
template
void transform_watershed<int32_t >(const int32_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                         size_t*   const arrayLbl,
                                         size_t          nDim,
                                         int32_t         threshold,
                                         int32_t         seedDistMin,
                                         eWaterDir       water,
                                         std::vector< std::set<size_t> >& neighborList);
template
void transform_watershed<uint32_t>(const uint32_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                         size_t*   const arrayLbl,
                                         size_t          nDim,
                                         uint32_t        threshold,
                                         uint32_t        seedDistMin,
                                         eWaterDir       water,
                                         std::vector< std::set<size_t> >& neighborList);
template
void transform_watershed<int64_t >(const int64_t*  const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                         size_t*   const arrayLbl,
                                         size_t          nDim,
                                         int64_t         threshold,
                                         int64_t         seedDistMin,
                                         eWaterDir       water,
                                         std::vector< std::set<size_t> >& neighborList);
template
void transform_watershed<uint64_t>(const uint64_t* const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                         size_t*   const arrayLbl,
                                         size_t          nDim,
                                         uint64_t        threshold,
                                         uint64_t        seedDistMin,
                                         eWaterDir       water,
                                         std::vector< std::set<size_t> >& neighborList);
template
void transform_watershed<float   >(const float*    const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                         size_t*   const arrayLbl,
                                         size_t          nDim,
                                         float           threshold,
                                         float           seedDistMin,
                                         eWaterDir       water,
                                         std::vector< std::set<size_t> >& neighborList);
template
void transform_watershed<double  >(const double*   const arraySrc, size_t nRow, size_t nCol, size_t nSec,
                                         size_t*   const arrayLbl,
                                         size_t          nDim,
                                         double          threshold,
                                         double          seedDistMin,
                                         eWaterDir       water,
                                         std::vector< std::set<size_t> >& neighborList);

} // namespace gem

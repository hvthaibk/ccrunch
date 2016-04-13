/***********************************************************************
 *  File:       cPolygon.cpp
 *
 *  Purpose:    Implementation of a polygon class
 *
 *  Author:     Lucia Martin Reixach, Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cPolygon.hpp"
#include <utility>

namespace gem {

template <typename T>
cPolygon<T>::cPolygon(const cPolygon& other)
    : cSetPoint3<T>(other)
{}

template <typename T>
bool cPolygon<T>::isPolygon(void) const
{
    bool    returnVal = (cSetPoint3<T>::at(0) == cSetPoint3<T>::at(cSetPoint3<T>::size()-1));

    for (size_t i = 0; i < cSetPoint3<T>::size(); i++) {
        if (cSetPoint3<T>::at(i)[2] != 0) {
            returnVal = false;
        }
    }

    return returnVal;
}

template <typename T>
bool cPolygon<T>::isSimple(void) const
{
    cSegment3<T>    s1, s2;

    for (size_t i = 0; i < cSetPoint3<T>::size()-1; i++) {
        s1[0] = cSetPoint3<T>::at(i);
        s1[1] = cSetPoint3<T>::at(i+1);

        for (size_t j = i+1; j < cSetPoint3<T>::size(); j++) {
            s2[0] = cSetPoint3<T>::at(j);
            s2[1] = cSetPoint3<T>::at((j+1) % cSetPoint3<T>::size());
            if (s1.isIntersect(s2)) {
                return false;
            }
        }
    }
    return true;
}

template <typename T1, class T2>
bool orderPair(const std::pair<T1, T2>& s1, const std::pair<T1, T2>& s2) { return (s1.first <= s2.first); }

template <typename T>
void cPolygon<T>::opFill(const cVector3<size_t>& imSize, cData3<T>& imData) const
{
    cSegment3<T> scanLine = cSegment3<T>(cVector3<T>(cPolygon<T>::getCoordMinX()-1, 0, 0), cVector3<T>(cPolygon<T>::getCoordMaxX()+1, 0, 0));
    cSegment3<T> *pSegments = cPolygon<T>::getSegments();
    std::vector<std::pair<T, char> > coordX;
    bool flag;

    for (size_t i = (size_t)round(cPolygon<T>::getCoordMaxY()); i > (size_t)round(cPolygon<T>::getCoordMinY()); i--)
    {
        scanLine[0][1] = (T)i;
        scanLine[1][1] = (T)i;
        flag = true;

        for (size_t j = 0; j < cPolygon<T>::size(); j++)
        {

            if (scanLine.isIntersect(pSegments[j]))
            {
                coordX.push_back(std::make_pair((size_t)round(scanLine.getIntersectPoint(pSegments[j])[0]), '-'));
            }
            else if (scanLine.isTouched(pSegments[j]))
            {
                size_t coordYInterception = (size_t)round(scanLine.getIntersectPoint(pSegments[j])[1]);
                if (coordYInterception <= pSegments[j][0][1] && coordYInterception <= pSegments[j][1][1])
                    coordX.push_back(std::make_pair((size_t)round(scanLine.getIntersectPoint(pSegments[j])[0]), 'm'));
                else
                    coordX.push_back(std::make_pair((size_t)round(scanLine.getIntersectPoint(pSegments[j])[0]), 'M'));
            }
        }

        std::sort(coordX.begin(), coordX.end(), orderPair<T,char>);

        if (coordX.size() > 0)
        {
            for (size_t j = 0; j < coordX.size() - 1; j++)
            {
                if (coordX[j].first == coordX[j+1].first && coordX[j].second != coordX[j+1].second)
                    coordX.erase(coordX.begin()+j);
            }

              for (size_t j = 0; j < coordX.size() - 1; j++)
            {
                if (flag)
                    for (size_t k = (size_t)round(coordX[j].first); flag ? k <= (size_t)round(coordX[j+1].first) : k < (size_t)round(coordX[j+1].first); k++)
                      if ((i-1) * imSize[1] + k < imData.getNelement())
                        imData[(i-1) * imSize[1] + k] = 0;
                flag = !flag;
            }

            coordX.clear();
        }
    }
}

template <typename T>
T cPolygon<T>::getCoordMinX(void) const
{
    T    minX;

    minX = cSetPoint3<T>::at(0)[0];
    for (size_t i = 1; i < cSetPoint3<T>::size(); i++)
        if (minX > cSetPoint3<T>::at(i)[0])
            minX = cSetPoint3<T>::at(i)[0];
    return minX;
}

template <typename T>
T cPolygon<T>::getCoordMaxX(void) const
{
    T    maxX;

    maxX = cSetPoint3<T>::at(0)[0];
    for (size_t i = 1; i < cSetPoint3<T>::size(); i++)
        if (maxX < cSetPoint3<T>::at(i)[0])
            maxX = cSetPoint3<T>::at(i)[0];
    return maxX;
}

template <typename T>
T cPolygon<T>::getCoordMinY(void) const
{
    T    minY;

    minY = cSetPoint3<T>::at(0)[1];
    for (size_t i = 1; i < cSetPoint3<T>::size(); i++)
        if (minY > cSetPoint3<T>::at(i)[1])
            minY = cSetPoint3<T>::at(i)[1];
    return minY;
}

template <typename T>
T cPolygon<T>::getCoordMaxY(void) const
{
    T    maxY;

    maxY = cSetPoint3<T>::at(0)[1];
    for (size_t i = 1; i < cSetPoint3<T>::size(); i++)
        if (maxY < cSetPoint3<T>::at(i)[1])
            maxY = cSetPoint3<T>::at(i)[1];
    return maxY;
}

template <typename T>
cSegment3<T> * cPolygon<T>::getSegments(void) const
{
    size_t numPoints = cPolygon<T>::size();
    cSegment3<T> *segments = new cSegment3<T>[numPoints];

    for (size_t i = 0;  i < numPoints; i++)
        segments[i] = cSegment3<T>(cSetPoint3<T>::at(i % numPoints), cSetPoint3<T>::at((i+1) % numPoints));

    return segments;
}


// instantiation
template class cPolygon<float >;
template class cPolygon<double>;

} // namespace gem

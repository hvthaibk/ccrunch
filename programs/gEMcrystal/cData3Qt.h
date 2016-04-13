/***********************************************************************
 *  File:       cData3Qt.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef CDATA3QT_H
#define CDATA3QT_H

#include "cData3.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QImage>
#pragma GCC diagnostic pop

using namespace gem;

template <class T>
class cData3Qt: public cData3<T>
{
public:
    cData3Qt(const cVector3<size_t>& size = cVector3<size_t>(0,0,0));
    cData3Qt(const cData3Qt& other);
    ~cData3Qt() {};
    void dataToQImage(QImage *, int);
};

#endif // CDATA3QT_H

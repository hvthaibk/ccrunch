/***********************************************************************
 *  File:       point.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef POINT_H
#define POINT_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QGraphicsItem>
#pragma GCC diagnostic pop

class Point: public QGraphicsItem
{
public:
    Point(int, int);
    void paint(QPainter *, const QStyleOptionGraphicsItem *, QWidget *);
    QRectF boundingRect() const {return QRectF(-7, -7, 15, 15); }
};

#endif // POINT_H

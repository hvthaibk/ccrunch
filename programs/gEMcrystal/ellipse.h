/***********************************************************************
 *  File:       ellipse.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef ELLIPSE_H
#define ELLIPSE_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QGraphicsEllipseItem>
#pragma GCC diagnostic pop

class Ellipse: public QGraphicsEllipseItem
{
public:
	Ellipse();
	Ellipse(QRectF);
	Ellipse(QRect);
	void paint(QPainter *, const QStyleOptionGraphicsItem *, QWidget *);
    QRectF boundingRect() const {return QRectF(-7, -7, 15, 15); }
};

#endif // ELLIPSE_H

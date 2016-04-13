/***********************************************************************
 *  File:       point.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "point.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QPainter>
#pragma GCC diagnostic pop

Point::Point(int coordX, int coordY): QGraphicsItem()
{
    setPos(coordX, coordY);
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemIgnoresTransformations);
}

void Point::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *)
{
    painter->setRenderHint(QPainter::Antialiasing);
    painter->setRenderHint(QPainter::Antialiasing);
    painter->setPen(QPen(QColor(0, 0, 200), 4));
	painter->setBrush(QBrush(QColor(0, 0, 200)));
	painter->drawEllipse(-2, -2, 4, 4);
	//painter->drawPoint(-1, -1);
}

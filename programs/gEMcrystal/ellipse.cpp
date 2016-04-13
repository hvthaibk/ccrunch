/***********************************************************************
 *  File:       ellipse.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "ellipse.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QPainter>
#pragma GCC diagnostic pop


Ellipse::Ellipse(): QGraphicsEllipseItem()
{
	setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemIgnoresTransformations);
	setPen(QPen(QColor(0, 0, 200), 4));
}

Ellipse::Ellipse(QRect rectE): QGraphicsEllipseItem(rectE)
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemIgnoresTransformations);
	setPen(QPen(QColor(0, 0, 200), 4));
}

Ellipse::Ellipse(QRectF rectFE): QGraphicsEllipseItem(rectFE)
{
    setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable | QGraphicsItem::ItemIgnoresTransformations);
	setPen(QPen(QColor(0, 0, 200), 4));
}

void Ellipse::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
	painter->setRenderHint(QPainter::Antialiasing);
    painter->setRenderHint(QPainter::Antialiasing);
    painter->setPen(QPen(QColor(0, 0, 200), 4));
	painter->setBrush(QBrush(QColor(0, 0, 200)));
}

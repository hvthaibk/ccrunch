/***********************************************************************
 *  File:       micrograph.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "micrograph.h"
#include <iostream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
//
#pragma GCC diagnostic pop

Micrograph::Micrograph() : QGraphicsScene()
{
}

void Micrograph::mousePressEvent(QGraphicsSceneMouseEvent *sEvent)
{
	if (sEvent->button() == Qt::LeftButton)
		emit clickCoord((int)gem::round(sEvent->scenePos().x()), (int)gem::round(sEvent->scenePos().y()));
}


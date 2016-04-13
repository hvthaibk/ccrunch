/***********************************************************************
 *  File:       micrograph.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef MICROGRAPH_H
#define MICROGRAPH_H

#include "macro.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#pragma GCC diagnostic pop


class Micrograph: public QGraphicsScene
{
    Q_OBJECT

public:
    Micrograph();

signals:
    void clickCoord(int, int);

protected:
	void mousePressEvent(QGraphicsSceneMouseEvent *);
};

#endif // SCENE_H

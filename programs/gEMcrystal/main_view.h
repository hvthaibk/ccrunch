/***********************************************************************
 *  File:       main_view.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef MAIN_VIEW_H
#define MAIN_VIEW_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QAction>
#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsSceneContextMenuEvent>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsView>
#include <QMenu>
#include <QKeySequence>
#include <QShortcut>
#include <QTimeLine>
#include <QVector>
#include <QWheelEvent>
#pragma GCC diagnostic pop

//#include "scene.h"

class MainView: public QGraphicsView
{
    Q_OBJECT

public:
//    View(QWidget *parent = 0);
    MainView(QGraphicsScene *vScene); //, QWidget *parent);

private:
    int _numScheduledScalings;
    
private slots:
	void scalingTime(qreal);
	void animFinished();
	inline void zoomIn() { scale(1.1, 1.1); }
	inline void zoomOut() { scale(0.89, 0.89); }

protected:
    void wheelEvent(QWheelEvent *);

};

#endif // MAIN_VIEW_H

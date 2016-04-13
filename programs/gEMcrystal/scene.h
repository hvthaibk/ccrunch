/***********************************************************************
 *  File:       scene.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef SCENE_H
#define SCENE_H

#include "point.h"
#include "macro.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QAction>
#include <QGraphicsEllipseItem>
#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsSceneContextMenuEvent>
#include <QGraphicsSceneMouseEvent>
#include <QMenu>
#include <QVector>
#pragma GCC diagnostic pop


class Scene: public QGraphicsScene
{
    Q_OBJECT

public:
    Scene();
    QVector<QGraphicsEllipseItem *> ellipseVector;
    QVector<QGraphicsItem *> polygonVector;
	int vSize;
	QPointF startPoint, endPoint;

private:
	bool flagEllipse;
	bool flagPolygon;
	bool flagActive;
	bool flagEllipsePainting;
  int  mainViewHeight;
  int  mainViewWidth;

public:
	inline void setEllipsePainting(bool b)  { flagEllipsePainting = b; }
	inline bool isEllipsePainting(void)     { return flagEllipsePainting; }
	inline bool isEllipse(void) { return flagEllipse; }
	inline bool isPolygon(void) { return flagPolygon; }
	inline bool isActive(void)  { return flagActive; }
  inline void setMainViewSize(int x, int y) { mainViewWidth = x; mainViewHeight = y; }
    
public slots:
	inline void setEllipse(bool b) { flagEllipse = b; }
	inline void setPolygon(bool b) { flagPolygon = b; }
	inline void setActive(bool b)  { flagActive = b; }
	
signals:
    void vectorNotEmpty(bool);

protected:
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *);
	void mouseMoveEvent(QGraphicsSceneMouseEvent *);
	void mouseReleaseEvent(QGraphicsSceneMouseEvent *);
	void contextMenuEvent(QGraphicsSceneContextMenuEvent *);
};

#endif // SCENE_H

/***********************************************************************
 *  File:       scene.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "scene.h"
#include <iostream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
//
#pragma GCC diagnostic pop

#define vSize 1 //vector size

QGraphicsEllipseItem *ellipseItem = new QGraphicsEllipseItem();

Scene::Scene() : QGraphicsScene()
{
	ellipseVector.reserve(vSize);
    polygonVector.reserve(vSize);
    addLine(0, 0, 0, 1, QPen(Qt::transparent, 1));
    ellipseItem->setPen(QPen(QColor(0, 0, 200), 4));
    emit vectorNotEmpty(false);
    setActive(false);
}

void Scene::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *sEvent)
{
	if (isActive())
	{
		if (isPolygon())
		{
			int x = (int) gem::round(sEvent->scenePos().x());
			int y = (int) gem::round(sEvent->scenePos().y());
      
      if (x < 0)
        x = 0;
      if (x >= mainViewWidth)
        x = mainViewWidth-1;
      if (y < 0)
        y = 0;
      if (y >= mainViewHeight)
        y = mainViewHeight-1;
        
			Point *point = dynamic_cast<Point*>(itemAt(x, y));

			if (point == 0 && sEvent->button() == Qt::LeftButton)
			{
				QGraphicsItem *item = new Point(x, y);
				addItem(item);
				polygonVector << item;
				emit vectorNotEmpty(true);
			}
		}
		else if (isEllipse())
		{
			if (sEvent->button() == Qt::LeftButton)
			{
				startPoint = endPoint = sEvent->scenePos();
				ellipseItem->setRect(QRectF(2*startPoint-endPoint, endPoint));
				addItem(ellipseItem);
				setEllipsePainting(true);
			}
		} 

		QGraphicsScene::mouseDoubleClickEvent(sEvent);
	}
}

void Scene::mouseMoveEvent(QGraphicsSceneMouseEvent *sEvent)
{
	if (isActive())
	{
		if (isEllipse() && isEllipsePainting())
		{
			endPoint = sEvent->scenePos();
			if ((startPoint != endPoint) && (sEvent->buttons() == Qt::LeftButton))
				{
					removeItem(ellipseItem);
					ellipseItem->setRect(QRectF(2*startPoint-endPoint, endPoint));
					addItem(ellipseItem);
				}
		}

		QGraphicsScene::mouseMoveEvent(sEvent);
	}		
}

void Scene::mouseReleaseEvent(QGraphicsSceneMouseEvent *sEvent)
{
	if (isActive())
	{
		if (isEllipse() && isEllipsePainting() && startPoint != endPoint)
			if (sEvent->button() == Qt::LeftButton)
			{
				QGraphicsEllipseItem *finishedEllipseItem = new QGraphicsEllipseItem(ellipseItem->rect());
				setEllipsePainting(false); 
				finishedEllipseItem->setPen(QPen(QColor(0, 0, 200), 4));
				finishedEllipseItem->setBrush(QBrush(QColor(0, 0, 200, 50)));
				finishedEllipseItem->setFlags(QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable);
				removeItem(ellipseItem);
				addItem(finishedEllipseItem);
				ellipseVector << finishedEllipseItem;
			}

		QGraphicsScene::mouseReleaseEvent(sEvent);
	}
}

void Scene::contextMenuEvent(QGraphicsSceneContextMenuEvent *sEvent)
{
	if (isActive())
	{
		int x = (int) gem::round(sEvent->scenePos().x());
		int y = (int) gem::round(sEvent->scenePos().y());
		Point *point = dynamic_cast<Point*>(itemAt(x, y));
		QGraphicsEllipseItem *ellipse = dynamic_cast<QGraphicsEllipseItem*>(itemAt(x, y));
		if (point == 0 && ellipse == 0)
			return;
		else if (point != 0)
		{
			if (point->flags() & QGraphicsItem::ItemIsSelectable)
			{
				QMenu menu;
				QAction *deleteAction = menu.addAction("Delete Point");
				if (menu.exec(sEvent->screenPos()) == deleteAction)
				{
					removeItem(point);
					polygonVector.remove(polygonVector.indexOf(point));
					delete point;
					if (polygonVector.empty())
						emit vectorNotEmpty(false);
				}
			}
		}
		else if (ellipse != 0)
		{
			if (ellipse->flags() & QGraphicsItem::ItemIsSelectable)
			{
				QMenu menu;
				QAction *deleteAction = menu.addAction("Delete Ellipse");
				if (menu.exec(sEvent->screenPos()) == deleteAction)
				{
					removeItem(ellipse);
					ellipseVector.remove(ellipseVector.indexOf(ellipse));
					delete ellipse;
				}
			}
		}
	}
}






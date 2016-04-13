/***********************************************************************
 *  File:       main_view.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "main_view.h"
#include "cMRC.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
//
#pragma GCC diagnostic pop

MainView::MainView(QGraphicsScene *vScene) : QGraphicsView(vScene)
{
	_numScheduledScalings = 0;

	QShortcut *zmIn = new QShortcut(QKeySequence("Ctrl++"), this);
	connect(zmIn, SIGNAL(activated()), this, SLOT(zoomIn()));
	QShortcut *zmOut = new QShortcut(QKeySequence("Ctrl+-"), this);
	connect(zmOut, SIGNAL(activated()), this, SLOT(zoomOut()));
}

void MainView::wheelEvent(QWheelEvent *sEvent)
{
	if (sEvent->buttons() == Qt::MiddleButton)
	{
		int numDegrees = sEvent->delta() / 8;
		int numSteps = numDegrees / 15;

		_numScheduledScalings += numSteps;

		if (_numScheduledScalings * numSteps < 0)  // if user moved the wheel in another direction, we reset previously scheduled scalings
			_numScheduledScalings = numSteps;

		QTimeLine *anim = new QTimeLine(350, this);
		anim->setUpdateInterval(20);

		connect(anim, SIGNAL(valueChanged(qreal)), this, SLOT(scalingTime(qreal)));
		connect(anim, SIGNAL(finished()), this, SLOT(animFinished()));
		anim->start();
	}
    QGraphicsView::wheelEvent(sEvent);
}

void MainView::scalingTime(qreal /* x */)
{
	qreal factor = 1.0 + qreal(_numScheduledScalings) / 300.0;
    scale(factor, factor);
}

void MainView::animFinished()
{
	if (_numScheduledScalings > 0)
		_numScheduledScalings--;
	else
		_numScheduledScalings++;
	sender()->~QObject();
}

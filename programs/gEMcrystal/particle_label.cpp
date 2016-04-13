/***********************************************************************
 *  File:       particle_label.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "particle_label.h"
#include <iostream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
//
#pragma GCC diagnostic pop

ParticleLabel::ParticleLabel() : QLabel()
{
	deleted = 0;
	selected = false;
}

void ParticleLabel::mousePressEvent(QMouseEvent *sEvent)
{
	if (sEvent->button() == Qt::LeftButton)
		select();
}

void ParticleLabel::select()
{
	emit currentPikIndex(index);
	emit drawRect();
	if (selected == false)
		color();
}

void ParticleLabel::color()
{
	if (selected == false)
	{	
		QPixmap *pix = new QPixmap();
		*pix = *pixmap(); 
		QPainter painter(pix);
		painter.setBrush(QBrush(QColor(0, 0, 200, 30)));
		painter.drawRect(0, 0, size, size);
		painter.end();
		setPixmap(*pix);
	}
	selected = true;
}

void ParticleLabel::clearColor()
{
	if (selected == true)
	{
		setPixmap(originalPix);
		if (deleted == 1)
			right();
		else if (deleted == 2)
			wrong();
		selected = false;
	}
}
	
void ParticleLabel::right()
{
	QPixmap *pix = new QPixmap();
	*pix = *pixmap(); 
	QPainter painter(pix);
	painter.setPen(QPen(QColor(21, 158, 21), 0.04*size));
	painter.drawRect(0, 0, size, size);
	painter.end();
	setPixmap(*pix);
	if (deleted == 2)
		emit notDeletePik(index);
	
	delete pix;	
	deleted = 1;
}

void ParticleLabel::wrong()
{
	QPixmap *pix = new QPixmap();
	*pix = *pixmap(); 
	QPainter painter(pix);
	painter.setPen(QPen(QColor(215, 8, 0), 0.04*size));
	painter.drawRect(0, 0, size, size);
	painter.end();
	setPixmap(*pix);
	emit deletePik(index);

	delete pix;
	deleted = 2;
}

void ParticleLabel::backUpPix() //call AFTER having set pixmap
{
	originalPix = QPixmap();
	originalPix = *pixmap();
}

void ParticleLabel::setCorrVal(float cVal)
{
	QString tmp;
	int fontSize = (int) gem::round(0.07*size);
	tmp.setNum(cVal, 'f', 4);
	QPixmap *pix = new QPixmap();
	*pix = *pixmap(); 
	QPainter painter(pix);
	painter.setPen(QPen(QColor(0, 0, 200), 10));
	painter.setFont(QFont("Monospace", fontSize));
	painter.drawText((int) gem::round(fontSize/2.4), (int) gem::round(fontSize*1.8), tmp);
	painter.end();
	setPixmap(*pix);

	delete pix;
}


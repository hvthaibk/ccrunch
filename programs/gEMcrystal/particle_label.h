/***********************************************************************
 *  File:       particle_label.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef PARTICLE_LABEL_H
#define PARTICLE_LABEL_H

#include "macro.hpp"
#include "cData3Qt.h"
#include "cMRC.hpp"
#include "main_view.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QGroupBox>
#include <QImage>
#include <QLabel>
#include <QMessageBox>
#include <QGraphicsRectItem>
#include <QGraphicsScene>
#include <QPushButton>
#include <QVBoxLayout>
#pragma GCC diagnostic pop


class ParticleLabel: public QLabel
{
    Q_OBJECT

public:
    ParticleLabel();
	void select();
	void color();
	void clearColor();
	void right();
	void wrong();
	inline void setSize() { size = pixmap()->width(); }
	inline int getSize() { return size; }
	void backUpPix();
	void setCorrVal(float);

	void setMicrographPath(QString path) { micrographPath = path; }
	QString getMicrographPath() { return micrographPath; }
	void setCoordX(int cX) { coordX = cX; }
	void setCoordY(int cY) { coordY = cY; }
	int getCoordX() { return coordX; }
	int getCoordY() { return coordY; }
	void setIndex(int i) { index = i; }
	int getIndex() { return index; }

private:
    QString micrographPath;
	int coordX;
	int coordY;
	int index;
	QPixmap originalPix;
	int size;

	// flags
	int deleted; // 0 nothing, 1 no (right), 2 yes (wrong)
	bool selected;

signals:
	void currentPikIndex(int);
	void deletePik(int);
	void notDeletePik(int);
	void drawRect();

protected:
    void mousePressEvent(QMouseEvent *);
};

#endif // PARTICLE_LABEL_H

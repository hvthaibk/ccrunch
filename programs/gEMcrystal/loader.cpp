/***********************************************************************
 *  File:       loader.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "loader.h"
#include "macro.hpp"

using namespace gem;


Loader::Loader(): QWidget()
{
	createDisplay();
}

void Loader::createDisplay()
{
	layout = new QHBoxLayout();
	//layout->setSpacing(30);
	setLayout(layout);
		
	QLabel *parFileL = new QLabel();
	parFileL->setText(QString("Parameters file"));
	layout->addWidget(parFileL);
	
	parFileLE = new QLineEdit();
	layout->addWidget(parFileLE);
	
	QPushButton *parFilePB = new QPushButton();
	parFilePB->setIcon(style()->standardIcon(QStyle::SP_DirIcon));
	parFilePB->setToolTip("Choose File");
	connect(parFilePB, SIGNAL(clicked()), this, SLOT(chooseFile()));
	layout->addWidget(parFilePB);
	
	parameters = std::vector<QString> (14, QString(""));
}

void Loader::chooseFile()
{
	QString fileName = QFileDialog::getOpenFileName (this, tr("Get Parameters File"), QDir::currentPath());
	loadFile(&fileName);
}

void Loader::loadFile(QString *fileName)
{
	parFileLE->setText(*fileName);
	if (!fileName->isEmpty())
	{ 
		loadParameters(fileName);
		emit parametersAvailable();
	}
}

void Loader::loadParameters(QString *string)
{
	QFile file(*string);
	parameters = std::vector<QString> (15, QString(""));
 
    file.open(QIODevice::ReadOnly);
	QTextStream textS(&file);
    while (!textS.atEnd())
	{
        QString line = textS.readLine();
        if (line != QString(""))
			process_line(line);
    }
	file.close();
}

void Loader::process_line(QString line)
{
	if (line.startsWith("target"))
	{ line.remove(0, 7); parameters[0] = line; return; }
	
	if (line.startsWith("search"))
	{ line.remove(0, 7); parameters[1] = line; return; }
	
	if (line.startsWith("maskRef"))
	{ line.remove(0, 8); parameters[2] = line; return; }
	
	if (line.startsWith("maskSch"))
	{ line.remove(0, 8); parameters[3] = line; return; }
	
	if (line.startsWith("result"))
	{ line.remove(0, 7); parameters[4] = line; return; }
	
	if (line.startsWith("mode"))
	{ line.remove(0, 5); parameters[5] = line; return; }
	
	if (line.startsWith("micContrast"))
	{ line.remove(0, 12); parameters[6] = line; return; }
	
	if (line.startsWith("rotAngle"))
	{ line.remove(0, 9); parameters[7] = line; return; }
	
	if (line.startsWith("threshold"))
	{ line.remove(0, 10); parameters[8] = line; return; }

	if (line.startsWith("highThreshold")) // inverse thresholdHigh -> highThreshold because of conflict with plain threshold
	{ line.remove(0, 14); parameters[9] = line; return; }
	
	if (line.startsWith("maxPiks"))
	{ line.remove(0, 8); parameters[10] = line; return; }
	
	if (line.startsWith("sizePiks"))
	{ line.remove(0, 9); parameters[11] = line; return; }

	if (line.startsWith("distPiks"))
	{ line.remove(0, 9); parameters[12] = line; return; }
	
	if (line.startsWith("CPU"))
	{ line.remove(0, 4); parameters[13] = line; return; }

	if (line.startsWith("GPU"))
	{ line.remove(0, 4); parameters[14] = line; return; }
}


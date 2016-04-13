/***********************************************************************
 *  File:       loader.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef LOADER_H
#define LOADER_H

#include "macro.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QAction>
#include <QComboBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFile>
#include <QFileDialog>
#include <QGroupBox>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QMessageBox>
#include <QProcess>
#include <QPushButton>
#include <QSpinBox>
#include <QSplitter>
#include <QString>
#include <QTextEdit>
#include <QTextStream>
#pragma GCC diagnostic pop


class Loader: public QWidget
{
    Q_OBJECT

public:
    Loader();
    
    inline std::vector<QString> * getParameters() { return &parameters; }
    void loadFile(QString *);
  	
private:

	QHBoxLayout *layout;
	QLineEdit *parFileLE;
	std::vector<QString> parameters;

	void createDisplay();
	void loadParameters(QString *);
	void process_line(QString);
	
private slots:

	void chooseFile();
	
signals:	

	void parametersAvailable();
};

#endif // LOADER_H

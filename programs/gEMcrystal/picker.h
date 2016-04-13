/***********************************************************************
 *  File:       picker.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef PICKER_H
#define PICKER_H

#include "macro.hpp"
#include "output_text_edit.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QAction>
#include <QComboBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFile>
#include <QFileDialog>
#include <QFormLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QMessageBox>
#include <QProcess>
#include <QPushButton>
#include <QScrollArea>
#include <QScrollBar>
#include <QSpinBox>
#include <QSplitter>
#include <QString>
#include <QTextEdit>
#include <QTextStream>
#pragma GCC diagnostic pop


class Picker: public QGroupBox
{
    Q_OBJECT

public:
    Picker(QSplitter *split);

	void loadParameters(std::vector<QString> *);
  	
  	QGroupBox *outputBox;
  	    
	QLineEdit *targetLE;
  	QLineEdit *resultLE;

private:

	// general
	QVBoxLayout *layout;
	QProcess *gEMpicker;
	QProcess *helpPick;
		
	// settings
	QLineEdit *setFileLE;
	QLineEdit *searchLE;
	QLineEdit *maskRLE;
	QLineEdit *maskSLE;
	QComboBox *modeC;
	QComboBox *contrastC;
	QDoubleSpinBox *rotationS;
	QDoubleSpinBox *thresholdS;
	QDoubleSpinBox *thresholdHighS;
	QSpinBox *maxNumberS;
	QSpinBox *sizeS;
	QSpinBox *distS;
	QSpinBox *cpuS;
	QSpinBox *gpuS;
	
	// output
	QSplitter *outputP;
	//QGroupBox *outputBox;
	QTextEdit *outputTE;
	
	// displays		
	void createDisplay();
	void createSettingsDisplay();
	void createRunDisplay();
	void createOutputDisplay();
	
private slots:

	// picker
	void showPicker();
	void runPicker();
	void stopPicker();
	void helpPicker();
	
	// directories
	void chooseDirectoryT();
	void chooseDirectoryS();
	void chooseDirectoryMR();
	void chooseDirectoryMS();
	void chooseDirectoryR();
	
	// output
	void clearOutput();
	void saveOutput();
	void updateOutput();
	void helpOutput();
};

#endif // PICKER_H

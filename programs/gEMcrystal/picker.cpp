/***********************************************************************
 *  File:       picker.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "picker.h"
#include "macro.hpp"
#include "cMRC.hpp"

#define MAX_NUM_PART_PICKED 1000
#define SIZE_PART_PICKED 1000
#define CPU_MAX 1000
#define GPU_MAX 1000


using namespace gem;
using namespace std;

Picker::Picker(QSplitter *split = 0) : QGroupBox(split)
{
	outputP = split;
	createDisplay();
}

void Picker::createDisplay()
{
	setTitle("gEMpicker");
	setVisible(false);

  layout = new QVBoxLayout();
	layout->setSpacing(30);
	setLayout(layout);

	createSettingsDisplay();
	createRunDisplay();
	createOutputDisplay();
	//layout->addStretch();
}

void Picker::showPicker()
{
	setVisible(!isVisible());
}

void Picker::createSettingsDisplay()
{
	// settings box and layout
/*  QScrollArea *settingsScroll = new QScrollArea(this);
  settingsScroll->setWidgetResizable(true);
  //settingsScroll->show();
  QGroupBox *settingsBox = new QGroupBox(settingsScroll);
  settingsBox->setTitle("Parameters");
  layout->addWidget(settingsScroll);
  //settingsScroll->setWidget(settingsBox); */

  QScrollArea *settingsScroll = new QScrollArea();
  QGroupBox *settingsBox = new QGroupBox();
  settingsBox->setTitle("Parameters");
  settingsScroll->setWidgetResizable(true);
  settingsScroll->setWidget(settingsBox);


  layout->addWidget(settingsScroll);

/*  QGroupBox *settingsBox = new QGroupBox();
  settingsBox->setTitle("Parameters");
  layout->addWidget(settingsBox); */

	QVBoxLayout *settingsLayout = new QVBoxLayout();
  settingsLayout->setSpacing(30);
  settingsLayout->setSizeConstraint(QLayout::SetMinAndMaxSize);
  settingsBox->setLayout(settingsLayout);

	// setDirectory layout
	QVBoxLayout *setDirectoryLayout = new QVBoxLayout();
	setDirectoryLayout->setSpacing(5);
	settingsLayout->addLayout(setDirectoryLayout);

	QLabel *targetL = new QLabel();
	//targetL->setWordWrap(true);
	targetL->setText(QString("Target images"));
	setDirectoryLayout->addWidget(targetL);
	QHBoxLayout *targetLayout = new QHBoxLayout();
	targetLE = new QLineEdit();
	targetLayout->addWidget(targetLE);
	QPushButton *targetPB = new QPushButton();
	targetPB->setIcon(style()->standardIcon(QStyle::SP_DirIcon));
	targetPB->setToolTip("Choose Directory");
	connect(targetPB, SIGNAL(clicked()), this, SLOT(chooseDirectoryT()));
	targetLayout->addWidget(targetPB);
	setDirectoryLayout->addLayout(targetLayout);

	QLabel *searchL = new QLabel();
	//searchL->setWordWrap(true);
	searchL->setText(QString("Search images"));
	setDirectoryLayout->addWidget(searchL);
	QHBoxLayout *searchLayout = new QHBoxLayout();
	searchLE = new QLineEdit();
	searchLayout->addWidget(searchLE);
	QPushButton *searchPB = new QPushButton();
	searchPB->setIcon(style()->standardIcon(QStyle::SP_DirIcon));
	searchPB->setToolTip("Choose Directory");
	connect(searchPB, SIGNAL(clicked()), this, SLOT(chooseDirectoryS()));
	searchLayout->addWidget(searchPB);
	setDirectoryLayout->addLayout(searchLayout);

	QLabel *maskRL = new QLabel();
	//maskRL->setWordWrap(true);
	maskRL->setText(QString("Masks for target images"));
	setDirectoryLayout->addWidget(maskRL);
	QHBoxLayout *maskRLayout = new QHBoxLayout();
	maskRLE = new QLineEdit();
	maskRLayout->addWidget(maskRLE);
	QPushButton *maskRPB = new QPushButton();
	maskRPB->setIcon(style()->standardIcon(QStyle::SP_DirIcon));
	maskRPB->setToolTip("Choose Directory");
	connect(maskRPB, SIGNAL(clicked()), this, SLOT(chooseDirectoryMR()));
	maskRLayout->addWidget(maskRPB);
	setDirectoryLayout->addLayout(maskRLayout);

	QLabel *maskSL = new QLabel();
	//maskSL->setWordWrap(true);
	maskSL->setText(QString("Masks for search images"));
	setDirectoryLayout->addWidget(maskSL);
	QHBoxLayout *maskSLayout = new QHBoxLayout();
	maskSLE = new QLineEdit();
	maskSLayout->addWidget(maskSLE);
	QPushButton *maskSPB = new QPushButton();
	maskSPB->setIcon(style()->standardIcon(QStyle::SP_DirIcon));
	maskSPB->setToolTip("Choose Directory");
	connect(maskSPB, SIGNAL(clicked()), this, SLOT(chooseDirectoryMS()));
	maskSLayout->addWidget(maskSPB);
	setDirectoryLayout->addLayout(maskSLayout);

	QLabel *resultL = new QLabel();
	//resultL->setWordWrap(true);
	resultL->setText(QString("Result images"));
	setDirectoryLayout->addWidget(resultL);
	QHBoxLayout *resultLayout = new QHBoxLayout();
	resultLE = new QLineEdit();
	resultLayout->addWidget(resultLE);
	QPushButton *resultPB = new QPushButton();
	resultPB->setIcon(style()->standardIcon(QStyle::SP_DirIcon));
	resultPB->setToolTip("Choose Directory");
	connect(resultPB, SIGNAL(clicked()), this, SLOT(chooseDirectoryR()));
	resultLayout->addWidget(resultPB);
	setDirectoryLayout->addLayout(resultLayout);

/*
// block1 layout
	QBoxLayout *block1Layout = new QVBoxLayout();
	settingsLayout->addLayout(block1Layout);

  QHBoxLayout *modeHB = new QHBoxLayout();
	QLabel *modeL = new QLabel();
	modeL->setText(QString("Mode"));
  modeHB->addWidget(modeL);
	modeC = new QComboBox();
	modeC->addItem(QString("Compute Correlation"));
	modeC->addItem(QString("Perform Picking"));
  modeHB->addWidget(modeC);
  block1Layout->addLayout(modeHB);

  QHBoxLayout *contrastHB = new QHBoxLayout();
	QLabel *contrastL = new QLabel();
	contrastL->setText(QString("Micrograph contrast"));
  contrastHB->addWidget(contrastL);
	contrastC = new QComboBox();
	contrastC->addItem(QString("Negative Peaks"));
	contrastC->addItem(QString("Positive Peaks"));
  contrastHB->addWidget(contrastC);
	block1Layout->addLayout(contrastHB);
*/

////////////////////////////////////////////////////////////////////////

  // block1 layout
	QGridLayout *block1Layout = new QGridLayout();
	block1Layout->setColumnStretch(0, 0);
	block1Layout->setColumnStretch(1, 1);
	block1Layout->setVerticalSpacing(20);
	settingsLayout->addLayout(block1Layout);

	QLabel *modeL = new QLabel();
	//modeL->setWordWrap(true);
	modeL->setText(QString("Mode"));
	block1Layout->addWidget(modeL, 0, 0);
	modeC = new QComboBox();
	modeC->addItem(QString("Compute Correlation"));
	modeC->addItem(QString("Perform Picking - box only"));
	modeC->addItem(QString("Perform Picking - with particle images"));
	block1Layout->addWidget(modeC, 0, 1);

	QLabel *contrastL = new QLabel();
	//contrastL->setWordWrap(true);
	contrastL->setText(QString("Micrograph contrast"));
	block1Layout->addWidget(contrastL);
	contrastC = new QComboBox();
	contrastC->addItem(QString("Negative Peaks"));
	contrastC->addItem(QString("Positive Peaks"));
	block1Layout->addWidget(contrastC);

	// block2 layout
	QGridLayout *block2Layout = new QGridLayout();
	block2Layout->setColumnStretch(0, 1);
	block2Layout->setColumnStretch(1, 0);
	block2Layout->setVerticalSpacing(20);
	settingsLayout->addLayout(block2Layout);

	QLabel *rotationL = new QLabel();
	//rotationL->setWordWrap(true);
	rotationL->setText(QString("In-plane rotating angle (in degrees)"));
	block2Layout->addWidget(rotationL, 0, 0);
	rotationS = new QDoubleSpinBox();
	rotationS->setDecimals(1);
	rotationS->setRange(0, 360);
	rotationS->setSingleStep(0.1);
	block2Layout->addWidget(rotationS, 0, 1);

	QLabel *thresholdL = new QLabel();
	//thresholdL->setWordWrap(true);
	thresholdL->setText(QString("Threshold (low) value for picking"));
	block2Layout->addWidget(thresholdL, 1, 0);
	thresholdS = new QDoubleSpinBox();
	thresholdS->setDecimals(2);
	thresholdS->setRange(0.01, 1);
	thresholdS->setSingleStep(0.01);
	block2Layout->addWidget(thresholdS, 1, 1);

	QLabel *thresholdHighL = new QLabel();
	//thresholdL->setWordWrap(true);
	thresholdHighL->setText(QString("Threshold (high) value for picking"));
	block2Layout->addWidget(thresholdHighL, 2, 0);
	thresholdHighS = new QDoubleSpinBox();
	thresholdHighS->setDecimals(2);
	thresholdHighS->setRange(0.01, 1);
	thresholdHighS->setSingleStep(0.01);
	block2Layout->addWidget(thresholdHighS, 2, 1);

	QLabel *maxNumberL = new QLabel();
	//maxNumberL->setWordWrap(true);
	maxNumberL->setText(QString("Maximum number of particles picked")); // from each micrograph" (0 = no limit)"));
	block2Layout->addWidget(maxNumberL, 3, 0);
	maxNumberS = new QSpinBox();
	maxNumberS->setMaximum(MAX_NUM_PART_PICKED);
	block2Layout->addWidget(maxNumberS, 3, 1);

	QLabel *sizeL = new QLabel();
	//sizeL->setWordWrap(true);
	sizeL->setText(QString("Size of picked particles")); // (0 = use search image size)"));
	block2Layout->addWidget(sizeL, 4, 0);
	sizeS = new QSpinBox();
	sizeS->setMaximum(SIZE_PART_PICKED);
	block2Layout->addWidget(sizeS, 4, 1);

	QLabel *distL = new QLabel();
	//sizeL->setWordWrap(true);
	distL->setText(QString("Minimum distance between particles")); // (0 = use search image size)"));
	block2Layout->addWidget(distL, 5, 0);
	distS = new QSpinBox();
	distS->setMinimum(0);
	distS->setMaximum(424242);
	block2Layout->addWidget(distS, 5, 1);

	QLabel *cpuL = new QLabel();
	//cpuL->setWordWrap(true);
	cpuL->setText(QString("Number of CPU cores"));
	block2Layout->addWidget(cpuL, 6, 0);
	cpuS = new QSpinBox();
	cpuS->setMaximum(CPU_MAX);
	cpuS->setValue(1);
	block2Layout->addWidget(cpuS, 6, 1);

	QLabel *gpuL = new QLabel();
	//gpuL->setWordWrap(true);
	gpuL->setText(QString("Number of GPU cores"));
	block2Layout->addWidget(gpuL, 7, 0);
	gpuS = new QSpinBox();
	gpuS->setMaximum(GPU_MAX);
	block2Layout->addWidget(gpuS, 7, 1);


  settingsScroll->setMinimumWidth(settingsBox->minimumSizeHint().width() + 16);
  settingsBox->show();
	settingsScroll->show();
}

void Picker::chooseDirectoryT()
{
	QString fileName = QFileDialog::getExistingDirectory(this, tr("Get Directory (*.mrc)"), QDir::currentPath());
	targetLE->setText(fileName);
}

void Picker::chooseDirectoryS()
{
	QString fileName = QFileDialog::getExistingDirectory(this, tr("Get Directory (*.mrc)"), QDir::currentPath());
	searchLE->setText(fileName);
}

void Picker::chooseDirectoryMR()
{
	QString fileName = QFileDialog::getExistingDirectory(this, tr("Get Directory (*.tif)"), QDir::currentPath());
	maskRLE->setText(fileName);
}

void Picker::chooseDirectoryMS()
{
	QString fileName = QFileDialog::getExistingDirectory(this, tr("Get Directory (*.tif)"), QDir::currentPath());
	maskSLE->setText(fileName);
}

void Picker::chooseDirectoryR()
{
	QString fileName = QFileDialog::getExistingDirectory(this, tr("Get Directory"), QDir::currentPath());
	resultLE->setText(fileName);
}

void Picker::createRunDisplay()
{
	// run box
	QHBoxLayout *runLayout = new QHBoxLayout();
	layout->addLayout(runLayout);

	runLayout->addStretch();
	QPushButton *runPB = new QPushButton();
	runPB->setFixedSize(85, 85);
	runPB->setText("Run");
	connect(runPB, SIGNAL(clicked()), this, SLOT(runPicker()));

	runLayout->addWidget(runPB);
	runLayout->addStretch();
}

void Picker::createOutputDisplay()
{
	outputBox = new QGroupBox();
	outputBox->setFlat(true);
	outputBox->setVisible(false);
    outputP->insertWidget(1, outputBox);
    //QList<int> list;
    //list << 75 << 25;
    //outputP->setSizes(list);

    QVBoxLayout *outputLayout = new QVBoxLayout();
    outputBox->setLayout(outputLayout);

	QHBoxLayout *outputButtonsLayout = new QHBoxLayout();

	QPushButton *clearOutputPB = new QPushButton();
	clearOutputPB->setIcon(style()->standardIcon(QStyle::SP_FileIcon));
	clearOutputPB->setToolTip(QString("Clear Output"));
	connect(clearOutputPB, SIGNAL(clicked()), this, SLOT(clearOutput()));
	outputButtonsLayout->addWidget(clearOutputPB);

	QPushButton *saveOutputPB = new QPushButton();
	saveOutputPB->setIcon(style()->standardIcon(QStyle::SP_DriveFDIcon));
	saveOutputPB->setToolTip(QString("Save Output"));
	connect(saveOutputPB, SIGNAL(clicked()), this, SLOT(saveOutput()));
	outputButtonsLayout->addWidget(saveOutputPB);

	QPushButton *stopPickerPB = new QPushButton();
	stopPickerPB->setIcon(style()->standardIcon(QStyle::SP_BrowserStop));
	stopPickerPB->setToolTip(QString("Stop Picker"));
	connect(stopPickerPB, SIGNAL(clicked()), this, SLOT(stopPicker()));
	outputButtonsLayout->addWidget(stopPickerPB);

	QPushButton *helpPickerPB = new QPushButton();
	helpPickerPB->setIcon(style()->standardIcon(QStyle::SP_MessageBoxQuestion));
	helpPickerPB->setToolTip(QString("Help"));
	connect(helpPickerPB, SIGNAL(clicked()), this, SLOT(helpPicker()));
	outputButtonsLayout->addWidget(helpPickerPB);

	outputButtonsLayout->addStretch();
	outputLayout->addLayout(outputButtonsLayout);

    outputTE = new OutputTextEdit();
 	outputLayout->addWidget(outputTE);
}

void Picker::runPicker()
{
	QString dirTarget, dirSearch, dirMaskR, dirMaskS, dirResult;
	bool contrast;
	float rotationAngle, threshold, thresholdHigh;
	unsigned int mode, maxParticles, sizeParticles, distParticles, cpuCores, gpuCores;

	dirTarget	= targetLE->text();
	dirSearch	= searchLE->text();
	dirMaskR	= maskRLE->text();
	dirMaskS	= maskSLE->text();
	dirResult	= resultLE->text();

	mode		= modeC->currentIndex();
	contrast	= (contrastC->currentText() == QString("Negative Peaks")) ? 0 : 1;

	rotationAngle	= (float)rotationS->value();
	threshold		= (float)thresholdS->value();
	thresholdHigh	= (float)thresholdHighS->value();

	maxParticles	= (unsigned)maxNumberS->value();
	sizeParticles	= (unsigned)sizeS->value();
	distParticles	= (unsigned)distS->value();
	cpuCores		= (unsigned)cpuS->value();
	gpuCores		= (unsigned)gpuS->value();

	gEMpicker = new QProcess();
  QStringList executables;
	QString path;
	QStringList arguments;
	QString tmp;

  executables = QDir("./bin/").entryList(QStringList(QString("gEMpicker_*")), QDir::Files | QDir::Executable, QDir::Name);
  path = "./bin/" + (!executables.empty() ? executables.first() : "");

  if (path == "./bin/")
  {
    executables = QDir("./").entryList(QStringList(QString("gEMpicker*")), QDir::Files | QDir::Executable, QDir::Name);
    path = "./" + (!executables.empty() ? executables.first() : "");
  }

  if (path == "./bin/" || path == "./")
    QMessageBox::information(this, tr("gEMpicker"), tr("gEMpicker could not be found."));
	else
  {
    arguments << QString("--dirTgt=").append(dirTarget);
    arguments << QString("--dirSch=").append(dirSearch);
    if (dirMaskR != QString(""))
      arguments << QString("--dirMskRef=").append(dirMaskR);
    arguments << QString("--dirMskSch=").append(dirMaskS);
    arguments << QString("--dirRes=").append(dirResult);
    arguments << QString("--mode=").append(tmp.setNum(mode));
    arguments << QString("--angle2D=").append(tmp.setNum(rotationAngle));
    arguments << QString("--contrast=").append(tmp.setNum(contrast));
    arguments << QString("--thresh=").append(tmp.setNum(threshold));
    arguments << QString("--threshHigh=").append(tmp.setNum(thresholdHigh));
    arguments << QString("--nPickMax=").append(tmp.setNum(maxParticles));
    arguments << QString("--boxSize=").append(tmp.setNum(sizeParticles));
    arguments << QString("--boxDist=").append(tmp.setNum(distParticles));
    arguments << QString("--nCPU=").append(tmp.setNum(cpuCores));
    arguments << QString("--nGPU=").append(tmp.setNum(gpuCores));

    outputTE->insertPlainText("\n\nRunning " + path + " " + arguments.join(" ") + "\n");
    gEMpicker->start(path, arguments, QIODevice::ReadWrite | QIODevice::Text);

    connect(gEMpicker, SIGNAL(readyReadStandardOutput()), this, SLOT(updateOutput()));
    connect(gEMpicker, SIGNAL(readyReadStandardError()), this, SLOT(updateOutput()));
    connect(gEMpicker, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(updateOutput()));

    outputBox->setVisible(true);
  }
}

void Picker::stopPicker()
{
	int answer =  QMessageBox::question(outputBox, QString("Stop gEMpicker"),
										QString("Please confirm. Do you really want to stop gEMpicker?"),
										QMessageBox::Yes | QMessageBox::No, QMessageBox::No);
	if (answer == QMessageBox::Yes) gEMpicker->kill();
}

void Picker::updateOutput()
{
  	QString tmp = gEMpicker->readAllStandardOutput();

	if (tmp.contains("%", Qt::CaseInsensitive))
	{
		QString text = outputTE->toPlainText();
		text.chop((text.size()-1) - text.lastIndexOf(QChar('\n')));
		outputTE->clear();
		outputTE->setPlainText(text.append(tmp));
		outputTE->moveCursor(QTextCursor::End);
	}
	else
	{
		outputTE->insertPlainText(tmp);
		outputTE->insertPlainText(gEMpicker->readAllStandardError());
	}
}

void Picker::clearOutput()
{
	int answer =  QMessageBox::question(outputBox, QString("Clear Output"),
										QString("Please Confirm. Do you really want to clear the output?"),
										QMessageBox::Yes | QMessageBox::No, QMessageBox::No);
	if (answer == QMessageBox::Yes) outputTE->clear();
}

void Picker::saveOutput()
{
	ofstream file;
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Output"), QDir::currentPath(), tr("Logfile (*.log)"));

    file.open(fileName.toStdString().c_str(), std::ios_base::app);
	QString tmp = outputTE->toPlainText();
	file << tmp.toStdString();

	file.close();
}


void Picker::helpPicker()
{
	helpPick = new QProcess();
  QStringList executables;
	QString path;
	QStringList arguments;

 	executables = QDir("./bin/").entryList(QStringList(QString("gEMpicker_*")), QDir::Files | QDir::Executable, QDir::Name);
  path = "./bin/" + (!executables.empty() ? executables.first() : "");

  if (path == "./bin/")
  {
    executables = QDir("./").entryList(QStringList(QString("gEMpicker*")), QDir::Files | QDir::Executable, QDir::Name);
    path = "./" + (!executables.empty() ? executables.first() : "");
  }
  arguments << QString("--help");
	helpPick->start(path, arguments, QIODevice::ReadWrite | QIODevice::Text);

	connect(helpPick, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(helpOutput()));
}

void Picker::helpOutput()
{
	outputTE->insertPlainText(helpPick->readAllStandardOutput());
}

void Picker::loadParameters(std::vector<QString> *parameters)
{
	if ((*parameters)[0] != QString(""))
		targetLE->setText((*parameters)[0]);
	if ((*parameters)[1] != QString(""))
		searchLE->setText((*parameters)[1]);
	if ((*parameters)[2] != QString(""))
		maskRLE->setText((*parameters)[2]);
	if ((*parameters)[3] != QString(""))
		maskSLE->setText((*parameters)[3]);
	if ((*parameters)[4] != QString(""))
		resultLE->setText((*parameters)[4]);
	if ((*parameters)[5] != QString(""))
		modeC->setCurrentIndex((*parameters)[5].toInt());
	if ((*parameters)[6] != QString(""))
		contrastC->setCurrentIndex((*parameters)[6].toInt());
	if ((*parameters)[7] != QString(""))
		rotationS->setValue((*parameters)[7].toDouble());
	if ((*parameters)[8] != QString(""))
		thresholdS->setValue((*parameters)[8].toDouble());
	if ((*parameters)[9] != QString(""))
		thresholdHighS->setValue((*parameters)[9].toDouble());
	else
		thresholdHighS->setValue(1.0);
	if ((*parameters)[10] != QString(""))
		maxNumberS->setValue((*parameters)[10].toInt());
	if ((*parameters)[11] != QString(""))
		sizeS->setValue((*parameters)[11].toInt());
	if ((*parameters)[12] != QString(""))
		distS->setValue((*parameters)[12].toInt());
	if ((*parameters)[13] != QString(""))
		cpuS->setValue((*parameters)[13].toInt());
	if ((*parameters)[14] != QString(""))
		gpuS->setValue((*parameters)[14].toInt());
}


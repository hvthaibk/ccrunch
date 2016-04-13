/***********************************************************************
 *  File:       result.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "result.h"
#include "cMRC.hpp"

using namespace gem;

Result::Result(QSplitter *split = 0) : QGroupBox(split)
{
    splitter = split;
    createDisplay();
}

void Result::createDisplay()
{
    setTitle("Results Display");
    setVisible(false);

    layout = new QVBoxLayout();
    layout->setSpacing(30);
    setLayout(layout);

    createPiksDisplay();
    createSettingsDisplay();
    layout->addStretch();
}

void Result::showResult()
{
    setVisible(!isVisible());
}

void Result::createSettingsDisplay()
{
    QVBoxLayout *settingsLayout = new QVBoxLayout();
    settingsLayout->setSpacing(30);
    layout->addLayout(settingsLayout);

    // setDirectory layout
    QVBoxLayout *setDirectoryLayout = new QVBoxLayout();
    setDirectoryLayout->setSpacing(5);
    settingsLayout->addLayout(setDirectoryLayout);

    QLabel *dirMicrographsL = new QLabel();
    dirMicrographsL->setText(QString("Micrographs Directory"));
    setDirectoryLayout->addWidget(dirMicrographsL);
    QHBoxLayout *dirMicrographsLayout = new QHBoxLayout();
    dirMicrographsLE = new QLineEdit();
    dirMicrographsLayout->addWidget(dirMicrographsLE);
    QPushButton *dirMicrographsPB = new QPushButton();
    dirMicrographsPB->setIcon(style()->standardIcon(QStyle::SP_DirIcon));
    dirMicrographsPB->setToolTip(QString("Choose directory"));
    connect(dirMicrographsPB, SIGNAL(clicked()), this, SLOT(chooseDirectoryM()));
    dirMicrographsLayout->addWidget(dirMicrographsPB);
    setDirectoryLayout->addLayout(dirMicrographsLayout);
    QLabel *dirResultsL = new QLabel();
    dirResultsL->setText(QString("Results Directory"));
    setDirectoryLayout->addWidget(dirResultsL);
    QHBoxLayout *dirResultsLayout = new QHBoxLayout();
    dirResultsLE = new QLineEdit();
    dirResultsLayout->addWidget(dirResultsLE);
    QPushButton *dirResultsPB = new QPushButton();
    dirResultsPB->setIcon(style()->standardIcon(QStyle::SP_DirIcon));
    dirResultsPB->setToolTip(QString("Choose directory"));
    connect(dirResultsPB, SIGNAL(clicked()), this, SLOT(chooseDirectoryR()));
    dirResultsLayout->addWidget(dirResultsPB);
    setDirectoryLayout->addLayout(dirResultsLayout);

    QGridLayout *gridLayout = new QGridLayout();
    gridLayout->setColumnStretch(0, 1);
    gridLayout->setColumnStretch(1, 0);
    gridLayout->setVerticalSpacing(10);
    settingsLayout->addLayout(gridLayout);

    QLabel *quantL = new QLabel();
    quantL->setWordWrap(false);
    quantL->setText(QString("Number of images per row"));
    gridLayout->addWidget(quantL, 0, 0);
    quantS = new QSpinBox();
    quantS->setMinimum(1);
    quantS->setValue(6);
    quantS->setKeyboardTracking(false);
    connect(quantS, SIGNAL(valueChanged(int)), this, SLOT(updatePiksDisplayI(int)));
    gridLayout->addWidget(quantS, 0, 1);

    QLabel *minCorrVL = new QLabel();
    minCorrVL->setWordWrap(false);
    minCorrVL->setText(QString("Minimum correlation value"));
    gridLayout->addWidget(minCorrVL, 1, 0);
    minCorrVS = new QDoubleSpinBox();
    minCorrVS->setMinimum(-1);
    minCorrVS->setMaximum(1);
    minCorrVS->setSingleStep(0.05);
    minCorrVS->setValue(-1);
    minCorrVS->setKeyboardTracking(false);
    connect(minCorrVS, SIGNAL(valueChanged(double)), this, SLOT(updatePiksDisplayD(double)));
    gridLayout->addWidget(minCorrVS, 1, 1);

    QLabel *maxCorrVL = new QLabel();
    maxCorrVL->setWordWrap(false);
    maxCorrVL->setText(QString("Maximum correlation value"));
    gridLayout->addWidget(maxCorrVL, 2, 0);
    maxCorrVS = new QDoubleSpinBox();
    maxCorrVS->setMinimum(-1);
    maxCorrVS->setMaximum(1);
    maxCorrVS->setSingleStep(0.05);
    maxCorrVS->setValue(1);
    maxCorrVS->setKeyboardTracking(false);
    connect(maxCorrVS, SIGNAL(valueChanged(double)), this, SLOT(updatePiksDisplayD(double)));
    gridLayout->addWidget(maxCorrVS, 2, 1);

    QHBoxLayout *displayLayout = new QHBoxLayout();
    settingsLayout->addLayout(displayLayout);

    displayLayout->addStretch();
    QPushButton *displayPB = new QPushButton();
    displayPB->setFixedSize(85, 85);
    displayPB->setText("Display");
    connect(displayPB, SIGNAL(clicked()), this, SLOT(loadMicrographs()));
    displayLayout->addWidget(displayPB);
    displayLayout->addStretch();

    iteratorDisplayedFlag = false;
}

void Result::createMicrographIteratorDisplay()
{
    if (iteratorDisplayedFlag == false)
        layout->removeItem(layout->itemAt(layout->count()-1));
    iteratorDisplayedFlag = true;

    micrographIteratorBox = new QGroupBox();
    micrographIteratorBox->setTitle("Current Micrograph");
    layout->addWidget(micrographIteratorBox);

    QVBoxLayout *micrographIteratorLayout = new QVBoxLayout();
    micrographIteratorLayout->setSpacing(10);
    micrographIteratorBox->setLayout(micrographIteratorLayout);

    QHBoxLayout *micrographNameIteratorLayout = new QHBoxLayout();
    QPushButton *leftArrowPB = new QPushButton();
    leftArrowPB->setFixedSize(30, 30);
    leftArrowPB->setIcon(style()->standardIcon(QStyle::SP_ArrowLeft));
    leftArrowPB->setToolTip(QString("Go to the previous micrograph"));
    connect(leftArrowPB, SIGNAL(clicked()), this, SLOT(goLeft()));
    micrographNameIteratorLayout->addWidget(leftArrowPB);
    currentMicrographCB = new QComboBox();
    for (int i = 0; i < micrographFileNameList.size(); i++)
        currentMicrographCB->addItem(micrographFileNameList[i]);
    connect(currentMicrographCB, SIGNAL(currentIndexChanged(int)), this, SLOT(updateCurrentMicrograph(int)));
    micrographNameIteratorLayout->addWidget(currentMicrographCB);
    QPushButton *rightArrowPB = new QPushButton();
    rightArrowPB->setFixedSize(30, 30);
    rightArrowPB->setIcon(style()->standardIcon(QStyle::SP_ArrowRight));
    rightArrowPB->setToolTip(QString("Go to the next micrograph"));
    connect(rightArrowPB, SIGNAL(clicked()), this, SLOT(goRight()));
    micrographNameIteratorLayout->addWidget(rightArrowPB);
    micrographIteratorLayout->addLayout(micrographNameIteratorLayout);

    pikListTE = new QTextEdit();
    pikListTE->setReadOnly(true);
    micrographIteratorLayout->addWidget(pikListTE);

    changesLayout = new QHBoxLayout();
    applyChangesPB = new QPushButton();
    applyChangesPB->setFixedSize(135, 45);
    applyChangesPB->setText("Apply Changes");
    applyChangesPB->setToolTip(QString("Apply the changes made and update the files"));
    connect(applyChangesPB, SIGNAL(clicked()), this, SLOT(applyChanges()));
    changesLayout->addWidget(applyChangesPB);
    discardChangesPB = new QPushButton();
    discardChangesPB->setFixedSize(135, 45);
    discardChangesPB->setText("Discard Changes");
    discardChangesPB->setToolTip(QString("Discard the changes made"));
    connect(discardChangesPB, SIGNAL(clicked()), this, SLOT(discardChanges()));
    changesLayout->addWidget(discardChangesPB);
    micrographIteratorLayout->addLayout(changesLayout);

    micrographIteratorBox->show();
}

void Result::createPiksDisplay()
{
    particleScroll = new QScrollArea();
    QGroupBox *particlesBox = new QGroupBox();
    particlesBox->setFlat(true);
    QVBoxLayout *particlesLayout = new QVBoxLayout();
    piksLayout = new QGridLayout();
    piksLayout->setSpacing(1);
    splitter->insertWidget(0, particleScroll);
    QList<int> list;
    list << 100 << 0;
    splitter->setSizes(list);
    particleScroll->setWidgetResizable(true);
    particleScroll->setWidget(particlesBox);
    particlesBox->show();
    particleScroll->show();
    particlesBox->setLayout(particlesLayout);
    particlesLayout->addLayout(piksLayout);
    particlesLayout->addStretch();
}

void Result::loadMicrographs()
{
    dirMicrographs = dirMicrographsLE->text();
    dirResults = dirResultsLE->text();

    if (iteratorDisplayedFlag == true)
        delete micrographIteratorBox;

    if (dirMicrographs == QString("") || dirResults == QString(""))
    {
        if (iteratorDisplayedFlag == true) layout->addStretch();
        iteratorDisplayedFlag = false;
        QMessageBox::information(this, tr("Choose Directory"), tr("Choose a directory."));
        return;
    }
    micrographFileNameList = QDir(dirMicrographs).entryList(QStringList(QString("*.mrc")), QDir::Files, QDir::Name);
    if (micrographFileNameList.empty())
    {
        if (iteratorDisplayedFlag == true) layout->addStretch();
        iteratorDisplayedFlag = false;
        QMessageBox::information(this, tr("Choose Directory"), tr("No micrographs found. Choose another directory."));
        return;
    }
    windowOpenFlag = false;
    createMicrographIteratorDisplay();
    updateCurrentMicrograph(0);
    loadPiks(0);
}

void Result::loadPiks(int i)
{
    QString currentMicrographFileName;
    QStringList pikFileNameList;
    QString tmp;

    QLayoutItem *item;
    while ((item = piksLayout->takeAt(0)) != NULL)
    {
        delete item->widget();
        delete item;
    }
    pikParticleLabelVector.clear();
    pikInfoVector.clear();
    piksForDeletingIndexVector.clear();
    pikSelected = false;

    currentMicrographFileName = micrographFileNameList[i];
    currentMicrographFileName.chop(4);

    pikInfoVector.read((dirResults + "/pik_coord/" + currentMicrographFileName + ".txt").toStdString());

    tmp = dirResults + "/pik_ext/" + currentMicrographFileName;
    pikFileNameList = QDir(tmp).entryList(QStringList(QString("*.mrc")), QDir::Files, QDir::Name);
    micrographIteratorBox->show();

    piksCont = 0;

    cData3Qt<float> imageData;
    cMRC<float>().read((dirMicrographsLE->text() + "/" + currentMicrographFileName + ".mrc").toStdString(), imageData, NULL);
    imageData.autocontrast8bit(1);
    refNrow = imageData.getNrow();

    for (int j = 0; j < pikFileNameList.size(); j++)
    {
        if (!pikFileNameList[j].isEmpty())
        {
      QImage *image = new QImage();
      // cData3Qt<float> particleImageDisplayed;
      // read((tmp + "/" + pikFileNameList[j]).toStdString(), particleImageDisplayed, NULL);
      cData3Qt<float> particleData;
      particleData.opCrop(imageData,
                          cSize3(pikInfoVector[j]._boxSize,pikInfoVector[j]._boxSize,1),
                          cSize3(imageData.getNrow()-pikInfoVector[j]._y, pikInfoVector[j]._x, 0)
                                - cSize3(pikInfoVector[j]._boxSize/2,pikInfoVector[j]._boxSize/2,0));
      particleData.dataToQImage(image, 256);
      particleLabel = new ParticleLabel();
      particleLabel->setPixmap(QPixmap::fromImage(*image));
      particleLabel->setSize();
      particleLabel->setCorrVal(pikInfoVector[j]._corrVal);
      particleLabel->backUpPix();
      particleLabel->setScaledContents(true);
      particleLabel->setMaximumSize((int)(splitter->width()/(quantS->value() * 1.045)), (int)(splitter->width()/(quantS->value() * 1.045)));
      particleLabel->setMinimumSize((int)(splitter->width()/(quantS->value() * 1.045)), (int)(splitter->width()/(quantS->value() * 1.045)));
      if (pikInfoVector[j]._corrVal >= minCorrVS->value() && pikInfoVector[j]._corrVal <= maxCorrVS->value())
      {
        piksLayout->addWidget(particleLabel, (int)(piksCont / quantS->value()), piksCont % quantS->value());
        piksCont++;
      }
      particleLabel->setIndex(j);
      particleLabel->setMicrographPath(dirMicrographsLE->text() + "/" + currentMicrographFileName + ".mrc");
      particleLabel->setCoordX((int) pikInfoVector[j]._x);
      particleLabel->setCoordY((int) pikInfoVector[j]._y);
      connect(particleLabel, SIGNAL(deletePik(int)), this, SLOT(deletePik(int)));
      connect(particleLabel, SIGNAL(notDeletePik(int)), this, SLOT(notDeletePik(int)));
      connect(particleLabel, SIGNAL(currentPikIndex(int)), this, SLOT(setCurrentPikIndex(int)));
      connect(particleLabel, SIGNAL(currentPikIndex(int)), this, SLOT(clearSelection(int)));
      connect(particleLabel, SIGNAL(drawRect()), this, SLOT(drawRect()));
      pikParticleLabelVector.push_back(particleLabel);
        }
    }

    if (windowOpenFlag == true)
        delete micrographBox;

    windowOpenFlag = true;
    micrographBox = new QGroupBox();
    micrographLayout = new QVBoxLayout();

    QHBoxLayout *buttonsLayout = new QHBoxLayout();
    micrographLayout->addLayout(buttonsLayout);
    QPushButton *rightPB = new QPushButton();
    rightPB->setIcon(style()->standardIcon(QStyle::SP_DialogYesButton));
    connect(rightPB, SIGNAL(clicked()), this, SLOT(right()));
    buttonsLayout->addWidget(rightPB);
    QPushButton *wrongPB = new QPushButton();
    wrongPB->setIcon(style()->standardIcon(QStyle::SP_DialogNoButton));
    connect(wrongPB, SIGNAL(clicked()), this, SLOT(wrong()));
    buttonsLayout->addWidget(wrongPB);
    buttonsLayout->addStretch();
    QLabel *currentMicrographL = new QLabel();
    currentMicrographL->setText(micrographFileNameList[currentMicrographIndex]);
    buttonsLayout->addWidget(currentMicrographL);

    scene = new Micrograph();
    MainView *view = new MainView(scene);
    QGraphicsPixmapItem *itemP = new QGraphicsPixmapItem();
    QImage *image = new QImage();
    connect(scene, SIGNAL(clickCoord(int, int)), this, SLOT(findPik(int, int)));

    micrographBox->setLayout(micrographLayout);
    micrographLayout->addWidget(view);

    QDesktopWidget desktop;
    QRect mainScreenSize = desktop.availableGeometry(desktop.primaryScreen());
    micrographBox->resize((int) (0.91*(double)mainScreenSize.height()),
                          (int) (0.94*(double)mainScreenSize.height()));
    scene->addItem(itemP);

    imageData.dataToQImage(image, 256);
    itemP->setPixmap(QPixmap::fromImage(*image));
    view->scale(0.944*0.91*mainScreenSize.height()/(double)refNrow, 0.944*0.91*mainScreenSize.height()/(double)refNrow);

    drawAllRect();

    QStringList reducedPikFileNameList;
    for (int j = 0; j < pikFileNameList.size(); j++)
        if (pikInfoVector[j]._corrVal >= minCorrVS->value() && pikInfoVector[j]._corrVal <= maxCorrVS->value())
            reducedPikFileNameList.push_back(pikFileNameList[j]);
    pikListTE->setText(reducedPikFileNameList.join("\n"));
    reducedPikFileNameList.clear();

    micrographBox->setVisible(true);
}


void Result::updateCurrentMicrograph(int i)
{
    currentMicrographIndex = i;
    loadPiks(i);
}

void Result::loadParameters(std::vector<QString> *parameters)
{
    if ((*parameters)[0] != QString(""))
        dirMicrographsLE->setText((*parameters)[0]);
    if ((*parameters)[4] != QString(""))
        dirResultsLE->setText((*parameters)[4]);
}

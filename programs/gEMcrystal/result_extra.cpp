/***********************************************************************
 *  File:       result_extra.cpp
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

void Result::chooseDirectoryM()
{
    QString fileName = QFileDialog::getExistingDirectory(this, tr("Choose Directory (*.mrc)"), QDir::currentPath());
    dirMicrographsLE->setText(fileName);
}

void Result::chooseDirectoryR()
{
    QString fileName = QFileDialog::getExistingDirectory(this, tr("Choose Directory (*.mrc)"), QDir::currentPath());
    dirResultsLE->setText(fileName);
}

// arrows

void Result::goLeft()
{
    if (piksForDeletingIndexVector.size() > 0)
        if (confirmDiscardChanges() == QMessageBox::No)
            applyChanges();
    if (currentMicrographCB->currentIndex() > 0)
        currentMicrographCB->setCurrentIndex(currentMicrographCB->currentIndex()-1);
}

void Result::goRight()
{
    if (piksForDeletingIndexVector.size() > 0)
        if (confirmDiscardChanges() == QMessageBox::No)
            applyChanges();
    if (currentMicrographCB->currentIndex() < currentMicrographCB->count())
        currentMicrographCB->setCurrentIndex(currentMicrographCB->currentIndex()+1);
}

// changes

int Result::confirmApplyChanges()
{
    int answer =  QMessageBox::question(this, QString("Apply Changes"),
            QString("Please Confirm. Do you really want to apply the changes?"),
            QMessageBox::Yes | QMessageBox::No, QMessageBox::No);
    return answer;
}

void Result::applyChanges()
{
    if (confirmApplyChanges() == QMessageBox::Yes)
    {
        QString micrographFolderName = micrographFileNameList[currentMicrographIndex];
        micrographFolderName.chop(4);
        QString path = dirResults + "/pik_ext/" + micrographFolderName + "/";
        for (int i = 0; i < (int) piksForDeletingIndexVector.size(); i++)
        {
            QStringList fileNames;
            fileNames = QDir(path).entryList(QStringList(QString::fromStdString(pikInfoVector[piksForDeletingIndexVector[i]]._numStr) + "*"), QDir::Files, QDir::Name);
            for (int j = 0; j < fileNames.size(); j++)
                QFile::remove(path + fileNames[j]);
            pikInfoVector[piksForDeletingIndexVector[i]]._numStr = "0000";
        }
        piksForDeletingIndexVector.clear();

        int i = 0;
        while (i < (int) pikInfoVector.size())
        {
            if (pikInfoVector[i]._numStr == "0000")
                pikInfoVector.erase(pikInfoVector.begin()+i);
            else
                i++;
        }
        pikInfoVector.write(((dirResults + "/pik_coord/" + micrographFolderName + "__piks.txt").toStdString()));
        updateCurrentMicrograph(currentMicrographCB->currentIndex());
    }
}

int Result::confirmDiscardChanges()
{
    int answer =  QMessageBox::question(this, QString("Discard Changes"),
            QString("Do you want to discard the changes made?"),
            QMessageBox::Yes | QMessageBox::No, QMessageBox::No);
    return answer;
}

void Result::discardChanges()
{
    if (confirmDiscardChanges() == QMessageBox::Yes)
        updateCurrentMicrograph(currentMicrographCB->currentIndex());
}

// selection of piks

void Result::drawRect()
{
    int cX = pikParticleLabelVector[currentPikIndex]->getCoordX();
    int cY = pikParticleLabelVector[currentPikIndex]->getCoordY();
    int labelSize = pikParticleLabelVector[currentPikIndex]->getSize();

    if (scene->items().length() > piksCont + 1)
        scene->removeItem(rect);
    rect = new QGraphicsRectItem(cX-labelSize/2, ((int)refNrow-cY)-labelSize/2, labelSize, labelSize);
    rect->setPen(QPen(QColor(0, 0, 200), 8));
    scene->addItem(rect);
    // micrographBox->showNormal();
    micrographBox->raise();
}

void Result::drawAllRect()
{
    int cX, cY, labelSize;

  for (size_t i = 0; i < pikParticleLabelVector.size(); i++)
  {
    if (pikInfoVector[i]._corrVal >= minCorrVS->value() && pikInfoVector[i]._corrVal <= maxCorrVS->value())
    {
      cX = pikParticleLabelVector[i]->getCoordX();
      cY = pikParticleLabelVector[i]->getCoordY();
      labelSize = pikParticleLabelVector[i]->getSize();

      rect = new QGraphicsRectItem(cX-labelSize/2, ((int)refNrow-cY)-labelSize/2, labelSize, labelSize);
      rect->setPen(QPen(Qt::white, 8));
      scene->addItem(rect);
    }
  }

    micrographBox->showNormal();
    micrographBox->raise();
}

void Result::findPik(int coordX, int coordY)
{
    for (int i = 0; i < (int) pikInfoVector.size(); i++)
    {
        if (coordX >= (int) pikInfoVector[i]._x - 100 && coordX <= (int) pikInfoVector[i]._x + 100)
            if (coordY >= ((int)refNrow - (int) pikInfoVector[i]._y) - 100 && coordY <= ((int)refNrow - (int) pikInfoVector[i]._y) + 100)
        if (pikInfoVector[i]._corrVal >= minCorrVS->value() && pikInfoVector[i]._corrVal <= maxCorrVS->value())
          pikParticleLabelVector[i]->select();
    }
}

void Result::clearSelection(int i)
{
    for (int j = 0; j < (int) pikParticleLabelVector.size(); j++)
    {
    if (i != j)
        pikParticleLabelVector[j]->clearColor();
    }
}

// delete slots

void Result::deletePik(int i)
{
    //bool flag = false;
    for (int j = 0; j < (int) piksForDeletingIndexVector.size(); j++)
        if (piksForDeletingIndexVector[j] == i)
            return;
    piksForDeletingIndexVector.push_back(i);
}

void Result::notDeletePik(int i)
{
    for (int j = 0; j < (int) piksForDeletingIndexVector.size(); j++)
        if (piksForDeletingIndexVector[j] == i)
        {
            piksForDeletingIndexVector.erase(piksForDeletingIndexVector.begin()+j);
            return;
        }
}


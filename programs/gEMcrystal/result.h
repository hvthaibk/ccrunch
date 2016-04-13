/***********************************************************************
 *  File:       result.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef RESULT_H
#define RESULT_H

#include "macro.hpp"
#include "cData3Qt.h"
#include "particle_label.h"
#include "micrograph.h"
#include "cResultPicker.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QAction>
#include <QComboBox>
#include <QDesktopWidget>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QGroupBox>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QMessageBox>
#include <QProcess>
#include <QPushButton>
#include <QScrollArea>
#include <QSpinBox>
#include <QSplitter>
#include <QString>
#include <QTextEdit>
#include <QWidget>
#pragma GCC diagnostic pop


class Result: public QGroupBox
{
    Q_OBJECT

public:
    size_t  refNrow;
    
    Result(QSplitter *split);
    
    void loadParameters(std::vector<QString> *);
    
    void createDisplay();
    void createSettingsDisplay();
    void createPiksDisplay();
    void createMicrographIteratorDisplay();
    void showResult();
    
    QScrollArea *particleScroll;

private:

    // general
    QVBoxLayout *layout;
    QSplitter *splitter;
    QString dirMicrographs;
    QString dirResults;    

    // panel
    QLineEdit *dirMicrographsLE;
    QLineEdit *dirResultsLE;
    QSpinBox *quantS;
    QDoubleSpinBox *minCorrVS;
    QDoubleSpinBox *maxCorrVS;
    QGroupBox *micrographIteratorBox;
    bool iteratorDisplayedFlag;
    QComboBox *currentMicrographCB;
    QTextEdit *pikListTE;
    QHBoxLayout *changesLayout;

    // micrograph
    QGroupBox *micrographBox;
    bool windowOpenFlag;
    QVBoxLayout *micrographLayout;
    Micrograph *scene;
    QGraphicsRectItem *rect;
    QStringList micrographFileNameList;
    std::vector<QString> micrographVector;
    int currentMicrographIndex;

    //  piks
    QGridLayout *piksLayout;
    ParticleLabel *particleLabel;
    cResultPicker<float> pikInfoVector;
    std::vector<ParticleLabel *> pikParticleLabelVector;
    std::vector<int> piksForDeletingIndexVector;
    int currentPikIndex;
    bool pikSelected;
    int piksCont;

    //changes
    QPushButton *applyChangesPB;
    QPushButton *discardChangesPB;
    int confirmApplyChanges();
    int confirmDiscardChanges();
    
private slots:
    void loadMicrographs();
    void loadPiks(int i);
    void chooseDirectoryM();
    void chooseDirectoryR();
        
    // delete
    void deletePik(int);
    void notDeletePik(int);

    // display of piks
    inline void updatePiksDisplay (      ) { updateCurrentMicrograph(currentMicrographCB->currentIndex()); }
    inline void updatePiksDisplayI(int   ) { updateCurrentMicrograph(currentMicrographCB->currentIndex()); }
    inline void updatePiksDisplayD(double) { updateCurrentMicrograph(currentMicrographCB->currentIndex()); }
    
    // selection of piks
    void drawRect();
    void drawAllRect();
    void findPik(int, int);
    void clearSelection(int);
    void updateCurrentMicrograph(int);
    inline void setCurrentPikIndex(int i) { currentPikIndex = i; pikSelected = true; }
    
    // arrows
    void goLeft();
    void goRight();
    
    // changes
    void applyChanges();
    void discardChanges();
    inline void right() { if (pikSelected == true) pikParticleLabelVector[currentPikIndex]->right(); }
    inline void wrong() { if (pikSelected == true) pikParticleLabelVector[currentPikIndex]->wrong(); }
};

#endif // RESULT_H

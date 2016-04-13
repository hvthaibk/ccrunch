/***********************************************************************
 *  File:       main_window.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"

#include <QAction>
#include <QApplication>
#include <QButtonGroup>
#include <QComboBox>
#include <QDir>
#include <QDirIterator>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QGraphicsPixmapItem>
#include <QGraphicsPolygonItem>
#include <QGraphicsScene>
#include <QGroupBox>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QMainWindow>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QPixmap>
#include <QPolygonF>
#include <QPrinter>
#include <QPushButton>
#include <QRadioButton>
#include <QScrollArea>
#include <QSplitter>
#include <QStatusBar>
#include <QTextEdit>
#include <QToolBar>
#include <QToolButton>
#include <QVariant>
#include <QWidget>

#pragma GCC diagnostic pop

#include "main_view.h"
#include "scene.h"
#include "picker.h"
#include "cData3Qt.h"
#include "result.h"
#include "loader.h"

#include "macro.hpp"

QT_BEGIN_NAMESPACE

class MainWindow : public QMainWindow

{

    Q_OBJECT

public:

    MainWindow();
    Loader *loader;

private slots:

    void open();
    void print();
    void close();
    void zoomIn();
    void zoomOut();
    void normalSize();
    //void fitToWindow();
	void addPoints();
	void addEllipses();
	void savePoints();
	void clearPoints();
	void drawPolygon(QPolygonF);
	void saveToFile();
	void createMaskDisplay();
	void updateAddPointsState(bool b) { addPointsB->setEnabled(b); }
	void loadParameters();
  void about();
	//void createResultsDisplay();
	void onlyOneRightPanel(const char *);
	inline void showViewer(void)	{ onlyOneRightPanel("viewer"); }
	inline void showMask(void)		{ onlyOneRightPanel("mask"); }
	inline void showPicker(void)	{ onlyOneRightPanel("picker"); }
	inline void showResult(void)	{ onlyOneRightPanel("result"); }

private:

    void createActions();
    void createMenus();
	void createDisplay(QSplitter *);
    void updateActions(bool b);
    void scaleImage(double factor);

    QWidget *centralWidget;
	QVBoxLayout *mainLayout;
	QHBoxLayout *buttonsLayout;

	// menu
    QMenuBar *menuBar;
    QMenu *fileMenu;
    QMenu *viewMenu;
    QMenu *helpMenu;

    QToolBar * toolBar;

   // QStatusBar *statusBar;

    QSplitter *splitter;
    QSplitter *vSplitter;

    MainView *graphicsView;
    Scene *scene;
	QGraphicsPixmapItem *pItem;
    int zoomC;
    QSize tmpS;

	// box
	QGroupBox *mainBox;
    QGroupBox *groupsBox;
    QGroupBox *maskBox;
	QGroupBox *littleButtonsBox;
	Picker *pickerBox;
	Result *resultsBox;
	
	// text edit
	QTextEdit *textEdit;
	
	// push button
	QPushButton *closeFilePB;
	QPushButton *savePointsB;
	QPushButton *addPointsB;
	QPushButton *clearPointsB;
	QPushButton *savePolygonsMaskB;

	QVBoxLayout *generalLayout;

	QRadioButton *polygonRB;
	QRadioButton *ellipseRB;

	cData3Qt<float> imageDisplayed;
	cData3Qt<float> imageModified;

#ifndef QT_NO_PRINTER
    QPrinter printer;
#endif

	// action and actionTool
    QAction *openAct;
    QAction *printAct;
    QAction *closeAct;
    QAction *exitAct;
    QAction *zoomInAct;
    QAction *zoomOutAct;
    QAction *normalSizeAct;
    QAction *zoomInActT;
    QAction *zoomOutActT;
    QAction *createNewMaskActT;
    QAction *callPickerActT;
    QAction *displayResultsActT;
  
	// tool button
    QToolButton *zoomInB;
    QToolButton *zoomOutB;
	QToolButton *createNewMaskB;
	QToolButton *callPickerB;
	QToolButton *displayResultsB;

	// flag
	bool flagMainViewOpen;

};


QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H

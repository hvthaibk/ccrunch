/***********************************************************************
 *  File:       main_window.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "main_window.h"
#include "main_view.h"
#include "scene.h"
#include "picker.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QtGui>
#pragma GCC diagnostic pop


MainWindow::MainWindow()
{

	mainBox = new QGroupBox();
	mainBox->setFlat(true);
//	mainBox->setStyleSheet("QGroupBox[flat=true] {border:0;}\n QScrollArea {border:0;}"); //\n MainView {border:0;}
	mainLayout = new QVBoxLayout();
	mainBox->setLayout(mainLayout);
	buttonsLayout = new QHBoxLayout();
	mainLayout->addLayout(buttonsLayout);

  splitter = new QSplitter();
  splitter->setGeometry(QRect(10, 10, 600, 400));
  splitter->setOrientation(Qt::Horizontal);
	splitter->setStyleSheet("QGroupBox[flat=true] {border: 2px solid #d4d0c8;}\n QScrollArea {border:0;}"); //\n MainView {border:0;}

	vSplitter = new QSplitter();
	vSplitter->setOrientation(Qt::Vertical);
	splitter->insertWidget(0, mainBox);
	mainLayout->addWidget(vSplitter);
	
	scene = new Scene();
    graphicsView = new MainView(scene);
    //graphicsView->setBackgroundRole(QPalette::Dark);
    graphicsView->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    graphicsView->setAlignment(Qt::AlignCenter);//Qt::AlignLeft | Qt::AlignTop);
    //graphicsView->setGeometry(20, 20, 400, 400);
	vSplitter->insertWidget(0, graphicsView);
	
	scene->setEllipse(false);
    scene->setPolygon(false);
	
	//splitter->setStretchFactor(0, 1);
	QList<int> list;
    list << 75 << 25;
    splitter->setSizes(list);
    zoomC = 0;

    setCentralWidget(splitter);
	createDisplay(vSplitter);
    createActions();
    createMenus();
    flagMainViewOpen = false;	

    setWindowTitle(tr("gEMcrystal"));
    statusBar();
    resize(1200, 900);
    showMaximized();
    
    showViewer();
}


void MainWindow::createActions()
{
    openAct = new QAction(tr("&Open..."), this);
    openAct->setShortcut(tr("Ctrl+O"));
    connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

    printAct = new QAction(tr("&Print..."), this);
    printAct->setShortcut(tr("Ctrl+P"));
    printAct->setEnabled(false);
    connect(printAct, SIGNAL(triggered()), this, SLOT(print()));

    closeAct = new QAction(tr("&Close"), this);
    closeAct->setShortcut(tr("Ctrl+C"));
    connect(closeAct, SIGNAL(triggered()), this, SLOT(close()));

    exitAct = new QAction(tr("&Exit"), this);
    exitAct->setShortcut(tr("Ctrl+E"));
    connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));

    zoomInAct = new QAction(tr("Zoom &In (25%)"), this);
    zoomInAct->setShortcut(tr("Ctrl++"));
    zoomInAct->setEnabled(false);
    connect(zoomInAct, SIGNAL(triggered()), this, SLOT(zoomIn()));

    zoomOutAct = new QAction(tr("Zoom &Out (25%)"), this);
    zoomOutAct->setShortcut(tr("Ctrl+-"));
    zoomOutAct->setEnabled(false);
    connect(zoomOutAct, SIGNAL(triggered()), this, SLOT(zoomOut()));

    normalSizeAct = new QAction(tr("&Normal Size"), this);
    normalSizeAct->setShortcut(tr("Ctrl+S"));
    normalSizeAct->setEnabled(false);
    connect(normalSizeAct, SIGNAL(triggered()), this, SLOT(normalSize()));

    //fitToWindowAct = new QAction(tr("&Fit to Window"), this);
    //fitToWindowAct->setEnabled(false);
    //fitToWindowAct->setCheckable(true);
    //fitToWindowAct->setShortcut(tr("Ctrl+F"));
    //connect(fitToWindowAct, SIGNAL(triggered()), this, SLOT(fitToWindow()));

/*
    zoomInActT = new QAction(tr("+"), this);
    zoomInActT->setIcon(QIcon(QString("./icons/zoom-in.png")));
    zoomInActT->setEnabled(false);
    connect(zoomInActT, SIGNAL(triggered()), this, SLOT(zoomIn()));

    zoomOutActT = new QAction(tr("-"), this);
    zoomOutActT->setIcon(QIcon(QString("./icons/zoom-out.png")));
    zoomOutActT->setEnabled(false);
    connect(zoomOutActT, SIGNAL(triggered()), this, SLOT(zoomOut()));
*/
}

void MainWindow::createMenus()
{
	QPushButton *viewerPB = new QPushButton();
    viewerPB->setText("Viewer");
    viewerPB->setToolTip(QString("Open Viewer to display\n.mrc and .tif images"));
    viewerPB->setFixedSize(100, 40);
    connect(viewerPB, SIGNAL(clicked()), this, SLOT(showViewer()));
    buttonsLayout->addWidget(viewerPB);
    
    QPushButton *maskPB = new QPushButton();
    maskPB->setText("Mask");
    maskPB->setToolTip(QString("Open Mask to\ncreate new masks"));
    maskPB->setFixedSize(100, 40);
    connect(maskPB, SIGNAL(clicked()), this, SLOT(showMask()));
    buttonsLayout->addWidget(maskPB);
    
	QPushButton *pickerPB = new QPushButton();
    pickerPB->setText("Picker");
    pickerPB->setToolTip(QString("Open Picker to\ncall gEMPicker"));
    pickerPB->setFixedSize(100, 40);
    connect(pickerPB, SIGNAL(clicked()), this, SLOT(showPicker()));
    buttonsLayout->addWidget(pickerPB);
    
    QPushButton *resultPB = new QPushButton();
    resultPB->setText("Result");
    resultPB->setToolTip(QString("Open Result to display\nthe results from gEMpicker"));
    resultPB->setFixedSize(100, 40);
    connect(resultPB, SIGNAL(clicked()), this, SLOT(showResult())); 
    buttonsLayout->addWidget(resultPB);

	buttonsLayout->addStretch();  
    
	loader = new Loader();
	buttonsLayout->addWidget(loader);
	connect(loader, SIGNAL(parametersAvailable()), this, SLOT(loadParameters()));
  
  QPushButton *aboutPB = new QPushButton();
  aboutPB->setIcon(style()->standardIcon(QStyle::SP_MessageBoxQuestion));
  connect(aboutPB, SIGNAL(clicked()), this, SLOT(about()));
  buttonsLayout->addWidget(aboutPB);
	
	littleButtonsBox = new QGroupBox;
	littleButtonsBox->setFlat(true);
	QHBoxLayout *littleButtonsLayout = new QHBoxLayout();
	littleButtonsBox->setLayout(littleButtonsLayout);
    
    QPushButton *openFilePB = new QPushButton();
    openFilePB->setIcon(style()->standardIcon(QStyle::SP_FileDialogContentsView));
    openFilePB->setToolTip(QString("Open an image"));
    connect(openFilePB, SIGNAL(clicked()), this, SLOT(open()));
	littleButtonsLayout->addWidget(openFilePB);
	
	closeFilePB = new QPushButton();
    closeFilePB->setIcon(style()->standardIcon(QStyle::SP_DialogCloseButton));
    closeFilePB->setToolTip(QString("Close the image"));
    connect(closeFilePB, SIGNAL(clicked()), this, SLOT(close()));
	littleButtonsLayout->addWidget(closeFilePB);
    
    buttonsLayout->addWidget(littleButtonsBox);
        
/*
    zoomInB = new QToolButton;
    zoomInB->setText("+");
    zoomInB->setGeometry(0,0,10,20);
    zoomInB->setDefaultAction(zoomInActT);

    zoomOutB = new QToolButton;
    zoomOutB->setGeometry(0,0,10,20);
    zoomOutB->setDefaultAction(zoomOutActT);
*/  

  //  menuBar = new QMenuBar(this);
  //  menuBar->addMenu(fileMenu);
  //  menuBar->addMenu(viewMenu);
  //  menuBar->addMenu(helpMenu);
  //  menuBar->setGeometry(QRect(0, 0, 800, 23));

   // statusBar = new QStatusBar();

    //setMenuBar(menuBar);
    //setStatusBar(statusBar);
}

void MainWindow::createDisplay(QSplitter *split)
{
	groupsBox = new QGroupBox();
    generalLayout = new QVBoxLayout();
    groupsBox->setLayout(generalLayout);
    groupsBox->setFlat(true);
    splitter->insertWidget(1, groupsBox);
    	
	createMaskDisplay();
	pickerBox = new Picker(split);
	generalLayout->addWidget(pickerBox);
	resultsBox = new Result(split);
	generalLayout->addWidget(resultsBox);
}


void MainWindow::updateActions(bool b)
{
    zoomInAct->setEnabled(b);
    zoomOutAct->setEnabled(b);
    normalSizeAct->setEnabled(b);
  //  zoomInActT->setEnabled(b);
  //  zoomOutActT->setEnabled(b);
}

void MainWindow::onlyOneRightPanel(const char *box)
{
	QList<int> list;
	list << 75 << 25;
		
	if (strcmp(box, "viewer") == 0)
	{
		groupsBox->setVisible(false);
		scene->setActive(false);
		littleButtonsBox->setVisible(true);
		vSplitter->widget(0)->setVisible(false);
		vSplitter->widget(1)->setVisible(true);
		if (vSplitter->count() > 1)
			vSplitter->setSizes(list);
			
	}
	else if (strcmp(box, "mask") == 0)
	{
		groupsBox->setVisible(true);
		maskBox->setVisible(true); 
		if (flagMainViewOpen)
			scene->setActive(true); 
		pickerBox->setVisible(false);
		resultsBox->setVisible(false);
		littleButtonsBox->setVisible(true);
		vSplitter->widget(0)->setVisible(false);
		vSplitter->widget(1)->setVisible(true);
		if (vSplitter->count() > 1)
			vSplitter->setSizes(list);
	}
	else if (strcmp(box, "picker") == 0)
	{
		groupsBox->setVisible(true);
		pickerBox->setVisible(true); //!pickerBox->isVisible());	
		maskBox->setVisible(false);
		scene->setActive(false);
		resultsBox->setVisible(false);
		littleButtonsBox->setVisible(false);
		vSplitter->widget(0)->setVisible(false);
		vSplitter->widget(1)->setVisible(false);
	}
	else if (strcmp(box, "result") == 0)
	{
		groupsBox->setVisible(true);
		resultsBox->setVisible(true); //!resultsBox->isVisible());
		maskBox->setVisible(false);
		scene->setActive(false);
		pickerBox->setVisible(false);
		littleButtonsBox->setVisible(false);
		vSplitter->widget(0)->setVisible(true);
		vSplitter->widget(1)->setVisible(false);
		if (vSplitter->count() > 1)
			vSplitter->setSizes(list);
	}
}


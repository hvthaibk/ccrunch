/***********************************************************************
 *  File:       functions.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
//#include <QtGui>
#pragma GCC diagnostic pop

#include "main_window.h"
#include "scene.h"

#include "cMRC.hpp"
#include "cTIFF.hpp"
#include "cPolygon.hpp"

using namespace gem;
using namespace std;

void MainWindow::open()
{
    if (flagMainViewOpen)
        close();

    QImage *image = new QImage();

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());

    if (!fileName.isEmpty())
    {
        if (fileName.endsWith(".tif", Qt::CaseInsensitive) || fileName.endsWith(".tiff", Qt::CaseInsensitive))
        {
            cTIFF<float>().read(fileName.toStdString(), imageDisplayed);
            imageDisplayed.dataToQImage(image, 256);
        }
        else if (fileName.endsWith(".mrc", Qt::CaseInsensitive))
        {
            cMRC<float>().read(fileName.toStdString(), imageDisplayed, NULL);
            imageDisplayed.dataToQImage(image, 256);
        }
        else
        {
            QMessageBox::information(this, tr("Image Viewer"), tr("Format not supported"));
            return;
        }

        if (image->isNull())
        {
            QMessageBox::information(this, tr("Image Viewer"), tr("Cannot load %1.").arg(fileName));
            return;
        }

        imageModified.memAlloc(imageDisplayed.getSize());
        imageModified.memSetVal(255);

        pItem = new QGraphicsPixmapItem();
        scene->addItem(pItem);

        QPixmap tmp = QPixmap::fromImage(*image);

        pItem->setPixmap(tmp);
        graphicsView->fitInView(pItem, Qt::KeepAspectRatio);

        normalSize();
        tmpS = tmp.size();

        printAct->setEnabled(true);
        //fitToWindowAct->setEnabled(true);
        updateActions(true);

        imageModified.memFree();
        imageModified.memAlloc(imageDisplayed.getSize());
        imageModified.memSetVal(255);

        flagMainViewOpen = true;
        statusBar()->showMessage(fileName);
        scene->setMainViewSize(tmp.width(), tmp.height());
    if (maskBox->isVisible())
            scene->setActive(true);
    connect(polygonRB, SIGNAL(toggled(bool)), scene, SLOT(setPolygon(bool)));
        connect(ellipseRB, SIGNAL(toggled(bool)), scene, SLOT(setEllipse(bool)));
        connect(scene, SIGNAL(vectorNotEmpty(bool)), this, SLOT(updateAddPointsState(bool)));
    }
}

void MainWindow::close()
{
    if (!flagMainViewOpen)
        return;

    int answer =  QMessageBox::question(this, QString("Close file"),
                                        QString("Please confirm. Do you really want to close this file?"),
                                        QMessageBox::Yes | QMessageBox::No, QMessageBox::No);
    if (answer == QMessageBox::Yes)
    {
        Scene *newScene = new Scene();
        graphicsView = dynamic_cast<MainView *>(graphicsView);
        graphicsView->setScene(newScene);
        //graphicsView->setBackgroundRole(QPalette::Dark);
        graphicsView->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
        graphicsView->setAlignment(Qt::AlignCenter);
        newScene->setEllipse(scene->isEllipse());
        newScene->setPolygon(scene->isPolygon());
        newScene->setActive(scene->isActive());
        delete scene;
        scene = newScene;
        splitter->setStretchFactor(0, 1);
        imageDisplayed.memFree();
        imageModified.memFree();
        clearPoints();
        flagMainViewOpen = false;
        scene->setActive(false);
        statusBar()->clearMessage();

        updateActions(false);
    }
    else
        return;
}


void MainWindow::print()
{

/*    Q_ASSERT(imageLabel->pixmap());

#ifndef QT_NO_PRINTER
    QPrintDialog dialog(&printer, this);
    if (dialog.exec()) {
        QPainter painter(&printer);
        QRect rect = painter.viewport();
        QSize size = pItem->pixmap()->size();
        size.scale(rect.size(), Qt::KeepAspectRatio);
        painter.setViewport(rect.x(), rect.y(), size.width(), size.height());
        painter.setWindow(imageLabel->pixmap()->rect());
        painter.drawPixmap(0, 0, *imageLabel->pixmap());
    }

#endif
*/
}

void MainWindow::zoomIn()
{
    zoomC++;
    scaleImage(1.20);
}

void MainWindow::zoomOut()
{
    zoomC--;
    scaleImage(1/1.20);
}

void MainWindow::normalSize()
{
    if (zoomC < 0)
        scaleImage(pow(1.20,-zoomC));
    else
        scaleImage(pow(1/1.20, zoomC));
    zoomC = 0;
}

/*
void MainWindow::fitToWindow()
{
    bool fitToWindow = fitToWindowAct->isChecked();
    graphicsView->fitInView(scene->itemsBoundingRect(),Qt::KeepAspectRatio);
    if (!fitToWindow) {
        graphicsView->resize(tmpS);
    }

    updateActions();
}
*/

void MainWindow::about()
{
  QMessageBox::about(this, "About", "<p></p><center><p><h2>gEMcrystal</h2></p></center>"
               "<p>The graphical tool for picking particles in</p>"
               "<p>cryo-electron micrographs using gEMpicker.</p>"
               "<small><center><p>&copy; Luc&iacute;a N. Mart&iacute;n Reixach</p></center>"
               "<center><p>Inria Nancy Grand Est, 2013</p></center></small>");

/* refuses to close if set manually
 QMessageBox *msg = new QMessageBox(this);

  msg->setStandardButtons(0);
  msg->setWindowTitle("About");
  msg->setText("<p></p><center><p><h2>gEMcrystal</h2></p></center>"
               "<p>The graphical tool for picking particles in</p>"
               "<p>cryo-electron micrographs using <i>gEMpicker</i>.</p>"
               "<small><center><p>&copy; Luc&iacute;a N. Mart&iacute;n Reixach</p></center>"
               "<center><p>INRIA Nancy Grand-Est, 2013</p></center></small>");
  msg->exec();
*/
}

void MainWindow::scaleImage(double factor)
{
    Q_ASSERT(pItem->pixmap());
    graphicsView->scale(factor, factor);

    graphicsView->scrollBarWidgets(Qt::AlignCenter);
}

void MainWindow::loadParameters()
{
    pickerBox->loadParameters(loader->getParameters());
    resultsBox->loadParameters(loader->getParameters());
}




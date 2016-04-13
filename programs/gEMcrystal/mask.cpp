/***********************************************************************
 *  File:       mask.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "main_window.h"
#include "scene.h"
#include "cMRC.hpp"
#include "cTIFF.hpp"
#include "cPolygon.hpp"

using namespace gem;
using namespace std;


void MainWindow::createMaskDisplay()
{
    maskBox = new QGroupBox();
    maskBox->setTitle("Mask");
    maskBox->setVisible(false);
    QButtonGroup *maskButtonGroup = new QButtonGroup();
    ellipseRB = new QRadioButton();
    polygonRB = new QRadioButton();
    QVBoxLayout *maskLayout = new QVBoxLayout();

    QGroupBox *ellipseBox = new QGroupBox();
    ellipseBox->setTitle("Ellipse");

    ellipseRB->setText("Create Ellipse");
    ellipseRB->setToolTip(QString("Click in the center and drag to create\na new one, drag to move an existing\none, right-click to delete it"));
    ellipseRB->setAutoExclusive(true);
    connect(ellipseRB, SIGNAL(toggled(bool)), scene, SLOT(setEllipse(bool)));
    maskButtonGroup->addButton(ellipseRB);

    QVBoxLayout *ellipseLayout = new QVBoxLayout();
    ellipseLayout->addWidget(ellipseRB);

    ellipseBox->setLayout(ellipseLayout);
    maskLayout->addWidget(ellipseBox);

    QGroupBox *polygonBox = new QGroupBox();
    polygonBox->setTitle("Polygon");

    addPointsB = new QPushButton(tr("Build Polygon"), polygonBox);
    addPointsB->setToolTip(QString("Use the vertices to create a simple polygon"));
    addPointsB->setEnabled(false);
    connect(addPointsB, SIGNAL(clicked()), this, SLOT(addPoints()));
    connect(scene, SIGNAL(vectorNotEmpty(bool)), this, SLOT(updateAddPointsState(bool)));
    savePointsB = new QPushButton(tr("Save Coordenates to .txt File"), polygonBox);
    savePointsB->setToolTip(QString("Create a .txt file with the coordenates from the display"));
    //savePointsB->setGeometry(QRect(20, 130, 80, 26));
    connect(savePointsB, SIGNAL(clicked()), this,  SLOT(savePoints()));
    clearPointsB = new QPushButton(tr("Clear Coordenates"), polygonBox);
    clearPointsB->setToolTip(QString("Remove the coordenates from the display"));
    //clearPointsB->setGeometry(QRect(20, 170, 80, 26));
    connect(clearPointsB, SIGNAL(clicked()), this,  SLOT(clearPoints()));
    savePolygonsMaskB = new QPushButton(tr("Save Mask to .tif File"), polygonBox);
    savePolygonsMaskB->setToolTip(QString("Create a .tif file with the new mask"));
    connect(savePolygonsMaskB, SIGNAL(clicked()), this, SLOT(saveToFile()));

    //groupBox->setStyleSheet("groupBox { border: 3px solid gray; border-radius: 40px; background: white; }");
    textEdit = new QTextEdit(polygonBox);
    textEdit->setReadOnly(true);
    textEdit->setMaximumWidth(130);
    //textEdit->setGeometry(QRect(10, 200, 150, 150));
    textEdit->setCurrentFont(QFont("Monospace"));
    textEdit->setFontPointSize(10.0);

    polygonRB->setText("Create Polygon");
    polygonRB->setToolTip(QString("Double click to create a new vertex, drag to move and existing one, right-click to delete it"));
    polygonRB->setAutoExclusive(true);
    connect(polygonRB, SIGNAL(toggled(bool)), scene, SLOT(setPolygon(bool)));
    maskButtonGroup->addButton(polygonRB);

    QVBoxLayout *polygonTextLayout = new QVBoxLayout();
    polygonTextLayout->addWidget(polygonRB);
    polygonTextLayout->addWidget(textEdit);

    QVBoxLayout *polygonButtonLayout = new QVBoxLayout();
    polygonButtonLayout->addWidget(addPointsB);
    polygonButtonLayout->addWidget(clearPointsB);
    polygonButtonLayout->addStretch();
    polygonButtonLayout->addWidget(savePointsB);
    polygonButtonLayout->addWidget(savePolygonsMaskB);

    QHBoxLayout *polygonLayout = new QHBoxLayout();
    polygonLayout->addLayout(polygonTextLayout);
    polygonLayout->addLayout(polygonButtonLayout);

    polygonBox->setLayout(polygonLayout);
    maskLayout->addWidget(polygonBox);

    maskBox->setLayout(maskLayout);

    generalLayout->addWidget(maskBox);
}

cPolygon<float> fromQPolygonFTocPolygon(QPolygonF qPoly)
{
    cPolygon<float> cPoly;
    cPoly.resize(qPoly.size());

    for (int i = 0; i < qPoly.size(); i++)
    {
        cPoly[i][0] = (float)qPoly[i].x();
        cPoly[i][1] = (float)qPoly[i].y();
    }

    return cPoly;
}

void MainWindow::addPoints()
{
    QString tmpSt;
    QPolygonF polygon;

    cPolygon<float> cPoly = cPolygon<float>();

    tmpSt.append(QString("%1%2\n").arg("x", 6, QChar(' ')).arg("y", 6, QChar(' ')));

    for (int i = 0; i < scene->polygonVector.size(); i++)
    {
        Point *tmpP = dynamic_cast<Point*>(scene->polygonVector[i]);
        tmpSt.append(QString("%1%2\n").arg((int)tmpP->pos().x(), 6, 10, QChar(' ')).arg((int)tmpP->pos().y(), 6, 10, QChar(' ')));
        tmpP->setFlags(QGraphicsItem::ItemIgnoresTransformations);
        polygon << scene->polygonVector[i]->pos();
    }

    cPoly = fromQPolygonFTocPolygon(polygon);

    if (cPoly.isSimple())
    {
        textEdit->append(tmpSt);
      tmpSt.clear();
        drawPolygon(polygon);

      cPoly.opFill(imageModified.getSize(), imageModified); //cVector3<float>(0, 0, 0),

        scene->polygonVector.clear();
        updateAddPointsState(false);
    }
    else
    {
        for (int i = 0; i < scene->polygonVector.size(); i++)
              scene->removeItem(scene->polygonVector[i]);
        QMessageBox::information(this, tr("Build Polygon"), tr("The polygon is not simple.\nCreate a new one."));
        scene->polygonVector.clear();
        updateAddPointsState(false);
    }
    return;
}

void opFillEllipse(QGraphicsEllipseItem *ellipse, const cVector3<float>& /* origin */, const cVector3<size_t>& imSize, cData3<float>& imData)
{
    // (x-c)^2 / (a-c)^2 + (y-b)^2 / (d-b)^2 <= 1

    size_t a = (size_t) gem::round(ellipse->rect().topLeft().x());
    size_t b = (size_t) gem::round((ellipse->rect().bottomRight().y() + ellipse->rect().topLeft().y()) / 2);
    size_t c = (size_t) gem::round((ellipse->rect().bottomRight().x() + ellipse->rect().topLeft().x()) / 2);
    size_t d = (size_t) gem::round(ellipse->rect().topLeft().y());



    for (size_t x = a; x <= ellipse->rect().bottomRight().x(); x++)
    {
        for (size_t y = d; y <= ellipse->rect().bottomRight().y(); y++)
        {
            if ((float((x-c)*(x-c)) / float((a-c)*(a-c))) + (float((y-b)*(y-b)) / float((d-b)*(d-b))) <= 1)
                imData[y * imSize[1] + x] = 0;
        }
    }
}

void MainWindow::addEllipses()
{
    for (int i = 0; i < scene->ellipseVector.size(); i++)
    {
        opFillEllipse(scene->ellipseVector[i], cVector3<float>(0, 0, 0), imageModified.getSize(), imageModified);
        scene->ellipseVector[i]->setFlag(QGraphicsItem::ItemIsMovable, false);
        scene->ellipseVector[i]->setFlag(QGraphicsItem::ItemIsSelectable, false);
    }
    scene->ellipseVector.clear();
}


void MainWindow::savePoints()
{
    ofstream file;

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Points"), QDir::currentPath(), tr("Plain Text (*.txt)"));

    file.open(fileName.toStdString().c_str(), std::ios_base::app);
    QString tmp = textEdit->toPlainText();

    file << tmp.toStdString();

    file.close();
    return;
}

void MainWindow::clearPoints()
{
    textEdit->clear();
    return;
}

void MainWindow::drawPolygon(QPolygonF polygon)
{
    QGraphicsPolygonItem *polygonItem = new QGraphicsPolygonItem();

    polygonItem->setBrush(QBrush(QColor(0, 0, 200, 50)));
    polygonItem->setPen(QPen(QColor(0, 0, 200), 2));
    polygonItem->setPolygon(polygon);
    scene->addItem(polygonItem);
    return;
}

void MainWindow::saveToFile()
{
    addEllipses();

    QStringList list = statusBar()->currentMessage().split(QChar('/'));
    QString micrographName = list.last();
    micrographName.chop(4);

    if (!(*loader->getParameters())[2].isEmpty())
        cTIFF<float>().write(imageModified, (*loader->getParameters())[2].toStdString() + "/" + micrographName.toStdString() + "__msk_ref.tif");
    else
    {
        list.removeLast();
        list.removeLast();
        QString path = list.join("/");
        if (!QDir(path).exists("msk_ref"))
            QDir(path).mkdir("msk_ref");
        cTIFF<float>().write(imageModified, path.toStdString() + "/msk_ref/" + micrographName.toStdString() + "__msk_ref.tif");
    }
}

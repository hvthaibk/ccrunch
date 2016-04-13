/***********************************************************************
 *  File:       main.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QApplication>
#pragma GCC diagnostic pop

#include "main_window.h"
#include "loader.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QStringList args = a.arguments();

    QString paramPath;
    QRegExp regParams("(-p|--parameters)");

    for (int i = 1; i < args.size(); ++i)
    {
        if (regParams.indexIn(args.at(i)) != -1 )
        {
            if (i+1 < args.size())
            {
                paramPath = args.at(i+1);
                i++;
            }
        }
        else
        {
            qDebug("Unknown arg: %s", qPrintable(args.at(i)));
        }
    }

    MainWindow w;
    w.loader->loadFile(&paramPath);
    w.show();

    return a.exec();
}

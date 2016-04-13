/***********************************************************************
 *  File:       cData3Qt.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "cData3Qt.h"

template <class T>
cData3Qt<T>::cData3Qt(const cVector3<size_t>& size) : cData3<T>(size)
{}

template <class T>
cData3Qt<T>::cData3Qt(const cData3Qt& other) : cData3<T>(other)
{}

template <class T>
void cData3Qt<T>::dataToQImage(QImage *image, int colors) //colors 256 or 1024
{
	int width, height, length;
	width = (int) cArray3<T>::getNcol();
	height = (int) cArray3<T>::getNrow();
	length = (int) cArray3<T>::getNelement();
	
    QImage tmpI(width, height, QImage::Format_Indexed8);

    cArray3<T>::normalize8bit();
	tmpI.setNumColors(colors);
    
    for(int i = 0; i < colors; i++)
    {
        QRgb myColor = qRgb(i, i, i);
        tmpI.setColor(i, myColor);
    }
    for (int i = 0; i < length; i++)
    {
		int x = i % width;
        tmpI.setPixel(x, (i-x)/width, (uint)(cArray3<T>::_data[(size_t) i]));
    }

    *image = tmpI;
}

template class cData3Qt<float >;
template class cData3Qt<double>;

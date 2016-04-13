/***********************************************************************
 *  File:       output_text_edit.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "output_text_edit.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
//
#pragma GCC diagnostic pop

OutputTextEdit::OutputTextEdit() : QTextEdit()
{
	setReadOnly(true);
	setContextMenuPolicy(Qt::CustomContextMenu);
 	
 	QColor c;
 	c.setRgb(46, 52, 54);
	QPalette p = QPalette();
	p.setColor(QPalette::Base, c);
	setPalette(p);
	
	QFont f = QFont("Monospace");
	f.setStyleHint(QFont::TypeWriter);
	setCurrentFont(f);
	setFontPointSize(9);
	setTextColor(Qt::white);
}

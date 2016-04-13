/***********************************************************************
 *  File:       output_text_edit.h
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#ifndef OUTPUT_TEXT_EDIT_H
#define OUTPUT_TEXT_EDIT_H

#include "macro.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <QTextEdit>
#pragma GCC diagnostic pop


class OutputTextEdit: public QTextEdit
{
    Q_OBJECT

public:
    OutputTextEdit();
			
protected:
    inline void mousePressEvent(QMouseEvent *) { return; }
    inline void mouseDoubleClickEvent(QMouseEvent *) { return; }
    inline void contextMenuEvent(QMouseEvent *) { return; }
    inline void mouseReleaseEvent(QMouseEvent *) { return; }
};

#endif // OUTPUT_TEXT_EDIT_H

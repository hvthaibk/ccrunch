/***********************************************************************
 *  File:       cString.hpp
 *
 *  Purpose:    Header file for a string class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CSTRING_HPP__
#define __GEM_CSTRING_HPP__

#include "macro.hpp"

#include <string>

namespace gem {

class cString : public std::string
{
private:

public:
     cString(void) : std::string("") {};
     cString(const std::string& str) : std::string(str) {};
    ~cString() {};

    // empty checking
    void    requireEmpty   (const std::string& funcName) const { require(length() == 0, funcName + ": cString object is not empty"); }
    void    requireNonEmpty(const std::string& funcName) const { require(length() >  0, funcName + ": cString object is empty");     }

    // operator overloading
    cString&    operator= (const cString&     other);
    cString&    operator= (const std::string& other);

    // uppercase / lowercase
    void    toUpper(size_t i) { at(i) = (char) std::toupper(at(i)); };
    void    toLower(size_t i) { at(i) = (char) std::tolower(at(i)); };
    void    toUpper(void);
    void    toLower(void);
    cString toUpper(void) const;
    cString toLower(void) const;

    // trim begin / end
    cString trimBegin(void) const;
    cString trimEnd  (void) const;
    cString trimBoth (void) const;

    // insert begin / end
    void    insertBegin(char filChar, size_t strLen);
    void    insertEnd  (char filChar, size_t strLen);

    // int2Str
    const cString& int2Str(unsigned int value, size_t strLen, char filChar = '0');
};

} // namespace gem

#endif

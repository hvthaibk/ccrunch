/***********************************************************************
 *  File:       cString.cpp
 *
 *  Purpose:    Implementation of a string class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cString.hpp"

namespace gem {

cString& cString::operator=(const cString& other)
{
    assign(other);

    return *this;
}

cString& cString::operator=(const std::string& other)
{
    assign(other);

    return *this;
}

void cString::toUpper(void)
{
    const std::string   funcName("void cString::toUpper(void)");

    requireNonEmpty(funcName);

    for (size_t i = 0; i < size(); i++) {
        toUpper(i);
    }
}

void cString::toLower(void)
{
    const std::string   funcName("void cString::toLower(void)");

    requireNonEmpty(funcName);

    for (size_t i = 0; i < size(); i++) {
        toLower(i);
    }
}

cString cString::toUpper(void) const
{
    const std::string   funcName("cString cString::toUpper(void) const");

    requireNonEmpty(funcName);

    cString     strUpper;

    strUpper.resize(size());

    for (size_t i = 0; i < size(); i++) {
        strUpper.toUpper(i);
    }

    return strUpper;
}

cString cString::toLower(void) const
{
    const std::string   funcName("cString cString::toLower(void) const");

    requireNonEmpty(funcName);

    cString     strLower;

    strLower.resize(size());

    for (size_t i = 0; i < size(); i++) {
        strLower.toLower(i);
    }

    return strLower;
}

cString cString::trimBegin(void) const
{
    const std::string   funcName("cString cString::trimBegin(void) const");

    requireNonEmpty(funcName);

    std::size_t     indBegin = find_first_not_of(" \a\b\f\n\r\t\v");

    return std::string(*this, indBegin);
}

cString cString::trimEnd(void) const
{
    const std::string   funcName("cString cString::trimEnd(void) const");

    requireNonEmpty(funcName);

    std::size_t     indEnd   = find_last_not_of (" \a\b\f\n\r\t\v");

    return std::string(*this, 0, indEnd - 1);
}

cString cString::trimBoth(void) const
{
    const std::string   funcName("cString cString::trimBoth(void) const");

    requireNonEmpty(funcName);

    std::size_t     indBegin = find_first_not_of(" \a\b\f\n\r\t\v");
    std::size_t     indEnd   = find_last_not_of (" \a\b\f\n\r\t\v");

    if(indBegin == std::string::npos) // no non-spaces
        return std::string("");

    return std::string(*this, indBegin, indEnd - indBegin + 1);
}

void cString::insertBegin(char filChar, size_t strLen)
{
    const std::string funcName("void cString::insertBegin(char filChar, size_t strLen)");

    require(length() <= strLen, funcName + ": requested length is too short");

    std::string     strTmp(strLen, filChar);

    strTmp.replace(strTmp.length()-length(), length(), *this);

    assign(strTmp);
}

void cString::insertEnd(char filChar, size_t strLen)
{
    const std::string funcName("void cString::insertEnd(char filChar, size_t strLen)");

    require(length() <= strLen, funcName + ": requested length is too short");

    std::string     strTmp(strLen, filChar);

    strTmp.replace(0, length(), *this);

    assign(strTmp);
}

const cString& cString::int2Str(unsigned int value, size_t strLen, char filChar)
{
    const std::string funcName("const cString cString::int2Str("
                                    "unsigned int value, size_t strLen, char filChar)");

    *this = num2str(value);
    insertBegin(filChar, strLen);

    return *this;
}

} // namespace gem

/***********************************************************************
 *  File:       cFileSystem.hpp
 *
 *  Purpose:    Header file for a filesystem class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CFILESYSTEM_HPP__
#define __GEM_CFILESYSTEM_HPP__

#include "cString.hpp"

#include <vector>

namespace gem {

class cFileSystem : public cString
{
private:

public:
     cFileSystem(void) : cString("") {};
     cFileSystem(const std::string& fileName) : cString(fileName) {};
    ~cFileSystem() {};

    // empty checking
    void    requireEmpty   (const std::string& funcName) const { require(length() == 0, funcName + ": cFileSystem object is not empty"); }
    void    requireNonEmpty(const std::string& funcName) const { require(length() >  0, funcName + ": cFileSystem object is empty");     }

    // operator overloading
    cFileSystem&    operator= (const std::string& other);
    cFileSystem&    operator= (const cFileSystem& other);

    const std::string   getFilePath(void) const;
    const std::string   getFileRoot(void) const;
    const std::string   getFileExt (void) const;
    void                getFileList(const std::string&        fileExt,
                                    std::vector<std::string>& fileList,
                                    bool                      fullPath = false) const;

    bool    isExist  (void) const;
    void    dirDelete(const std::string& hostname = "") const;
    void    dirCreate(const std::string& hostname = "") const;
};

} // namespace gem

#endif

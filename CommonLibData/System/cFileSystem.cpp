/***********************************************************************
 *  File:       cFileSystem.cpp
 *
 *  Purpose:    Implementation of a filesystem class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cFileSystem.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/filesystem.hpp"
#pragma GCC diagnostic pop

namespace gem {

#ifdef __unix__

const std::string   FILESYSTEM_PATH_SEP = "/";
const std::string   FILESYSTEM_EXT_SEP  = ".";

#elif __win32__

const std::string   FILESYSTEM_PATH_SEP = "\\";
const std::string   FILESYSTEM_EXT_SEP  = ".";

#endif

cFileSystem& cFileSystem::operator=(const std::string& other)
{
    assign(other);

    return *this;
}

cFileSystem& cFileSystem::operator=(const cFileSystem& other)
{
    assign(other);

    return *this;
}

const std::string cFileSystem::getFilePath(void) const
{
    const std::string   funcName("const std::string cFileSystem::getFilePath(void) const");

    requireNonEmpty(funcName);

    std::string     filePath;
    size_t          idx = rfind(FILESYSTEM_PATH_SEP);

    if (idx != npos) {
        filePath = substr(0,idx+1);
    }
    else {
        filePath = "";
    }

    return filePath;
}

const std::string cFileSystem::getFileRoot(void) const
{
    const std::string   funcName("const std::string cFileSystem::getFileRoot(void) const");

    requireNonEmpty(funcName);

    std::string     fileRoot;
    size_t          idx1 = rfind(FILESYSTEM_PATH_SEP);
    size_t          idx2 = rfind(FILESYSTEM_EXT_SEP);

    if (idx2 != npos) {
        if (idx1 < npos) {
            fileRoot = substr(idx1+1,idx2-idx1-1);
        }
        else {
            fileRoot = substr(0,idx2-idx1-1);
        }
    }
    else {
        fileRoot = "";
    }

    return fileRoot;
}

const std::string cFileSystem::getFileExt(void) const
{
    const std::string   funcName("const std::string cFileSystem::getFileExt(void) const");

    requireNonEmpty(funcName);

    std::string     fileExt;
    size_t          idx = rfind(FILESYSTEM_EXT_SEP);

    if (idx != npos) {
        fileExt = substr(idx);
    }
    else {
        fileExt = "";
    }

    return fileExt;
}

void cFileSystem::getFileList(const std::string&        fileExt,
                              std::vector<std::string>& fileList,
                              bool                      fullPath) const
{
    const std::string   funcName("void cFileSystem::getFileList("
                                    "const std::string& fileExt, "
                                    "std::vector<std::string>& fileList"
                                    "bool fullPath) const");

    cString                                 fileExtLower(fileExt);
    cString                                 fileExtCur;
    boost::filesystem::directory_iterator   iterBegin, iterEnd;

    fileExtLower.toLower();
    iterBegin = boost::filesystem::directory_iterator(*this);

    for (; iterBegin != iterEnd; ++iterBegin) {

        if (is_regular_file(*iterBegin)) {

            fileExtCur = boost::filesystem::extension(*iterBegin);
            fileExtCur.toLower();

            if (fileExtCur == fileExtLower) {
                if (fullPath) {
                    fileList.push_back(iterBegin->path().string());
                }
                else {
                    fileList.push_back(iterBegin->path().filename().string());
                }
            }
        }
    }

    std::sort(fileList.begin(), fileList.end());
}

bool cFileSystem::isExist(void) const
{
    const std::string   funcName("bool cFileSystem::isExist(void) const");

    requireNonEmpty(funcName);

    return boost::filesystem::exists(*this);
}

void cFileSystem::dirDelete(const std::string& hostname) const
{
    const std::string   funcName("void cFileSystem::dirDelete("
                                    "const std::string& hostname) const");

    requireNonEmpty(funcName);

    if (boost::filesystem::exists(*this)) {
        if (hostname.length()) {
            require(boost::filesystem::remove_all(*this), funcName +
                    ": cannot delete directory " + *this
                    + " in (" + hostname + ")");
        }
        else {
            require(boost::filesystem::remove_all(*this), funcName +
                    ": cannot delete directory " + *this);
        }
    }
    else {
#ifndef NDEBUG
        WARNING(funcName, "trying to delete a non-existing directory " + *this);
#endif
    }
}

void cFileSystem::dirCreate(const std::string& hostname) const
{
    const std::string   funcName("void cFileSystem::dirCreate("
                                    "const std::string& hostname) const");

    requireNonEmpty(funcName);

    std::string     dirName;

    if (at(length()-1) == '/') {
        dirName = substr(0,length()-1);
    }
    else {
        dirName = *this;
    }

    if (!boost::filesystem::exists(dirName)) {
        if (hostname.length()) {
            require(boost::filesystem::create_directories(dirName),
                    funcName + ": cannot create directory " + dirName
                    + " in (" + hostname + ")");
        }
        else {
            require(boost::filesystem::create_directories(dirName),
                    funcName + ": cannot create directory " + dirName);
        }
    }
    else {
#ifndef NDEBUG
        WARNING(funcName, "trying to create an existing directory " + dirName);
#endif
    }
}

} // namespace gem

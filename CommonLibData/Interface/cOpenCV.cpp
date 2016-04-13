/***********************************************************************
 *  File:       cOpenCV.cpp
 *
 *  Purpose:    Implementation of an OpenCV class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cOpenCV.hpp"

#ifdef __GEM_USE_OPENCV__

#include <cstdarg>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wconversion"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core.hpp"
#pragma GCC diagnostic pop

namespace gem {

template <typename T>
void cOpenCV<T>::putText(const cData3<T>& object, const std::string& fileName,
                         const std::string& textStr,
                         eColorName textColor,
                         T scale,
                         eNorm8bit norm8bit)
{
    const std::string   funcName("void cOpenCV::putText("
                                    "const cData3<T>& object, "
                                    "const std::string& fileName, "
                                    "const std::string& textStr, "
                                    "eColorName textColor, "
                                    "T scale, eNorm8bit norm8bit)");

    require(object.getAddrData()  != NULL, funcName + ": writing an empty object to an OpenCV file");
    require(object.getDimension() == 2,    funcName + ": writing a non-image object to an OpenCV file");
    require(textStr.length()      != 0,    funcName + ": putting an empty text string onto an OpenCV image");

    cv::Mat         cvData     ((int) object.getNrow(), (int) object.getNcol(), CV_8U  );
    cv::Mat         cvDataColor((int) object.getNrow(), (int) object.getNcol(), CV_8UC3);
    cData3<T>       dataSave;

    switch (norm8bit) {
        case NORM8BIT_TRUE:
            dataSave = object;
            dataSave.normalize8bit();

            #pragma omp parallel for
            for(size_t i = 0; i < object.getNelement(); i++) {
                 cvData.data[i] = (uint8_t) dataSave[i];
            }

            cvtColor(cvData, cvDataColor, cv::COLOR_GRAY2RGB);
            break;
        case NORM8BIT_FALSE:
            #pragma omp parallel for
            for(size_t i = 0; i < object.getNelement(); i++) {
                 cvData.data[i] = (uint8_t) object[i];
            }

            cvtColor(cvData, cvDataColor, cv::COLOR_GRAY2RGB);
            break;
        default:
            ERROR(funcName, "unsupported norm8bit mode");
    }

    switch (textColor) {
        case COLOR_WHITE:
            cv::putText(cvDataColor, textStr, cv::Point(0,(int)(30*scale)),
                        cv::FONT_HERSHEY_SIMPLEX, scale, cv::Scalar(255,255,255), (int)(2*scale), cv::LINE_AA);
            break;
        case COLOR_BLACK:
            cv::putText(cvDataColor, textStr, cv::Point(0,(int)(30*scale)),
                        cv::FONT_HERSHEY_SIMPLEX, scale, cv::Scalar(  0,  0,  0), (int)(2*scale), cv::LINE_AA);
            break;
        case COLOR_RED:
            cv::putText(cvDataColor, textStr, cv::Point(0,(int)(30*scale)),
                        cv::FONT_HERSHEY_SIMPLEX, scale, cv::Scalar(  0,  0,255), (int)(2*scale), cv::LINE_AA);
            break;
        case COLOR_GREEN:
            cv::putText(cvDataColor, textStr, cv::Point(0,(int)(30*scale)),
                        cv::FONT_HERSHEY_SIMPLEX, scale, cv::Scalar(  0,255,  0), (int)(2*scale), cv::LINE_AA);
            break;
        case COLOR_BLUE:
            cv::putText(cvDataColor, textStr, cv::Point(0,(int)(30*scale)),
                        cv::FONT_HERSHEY_SIMPLEX, scale, cv::Scalar(255,  0,  0), (int)(2*scale), cv::LINE_AA);
            break;
        case COLOR_CYAN:
            cv::putText(cvDataColor, textStr, cv::Point(0,(int)(30*scale)),
                        cv::FONT_HERSHEY_SIMPLEX, scale, cv::Scalar(255,255,  0), (int)(2*scale), cv::LINE_AA);
            break;
        case COLOR_MAGNETA:
            cv::putText(cvDataColor, textStr, cv::Point(0,(int)(30*scale)),
                        cv::FONT_HERSHEY_SIMPLEX, scale, cv::Scalar(255,  0,255), (int)(2*scale), cv::LINE_AA);
            break;
        case COLOR_YELLOW:
            cv::putText(cvDataColor, textStr, cv::Point(0,(int)(30*scale)),
                        cv::FONT_HERSHEY_SIMPLEX, scale, cv::Scalar(  0,255,255), (int)(2*scale), cv::LINE_AA);
            break;
        default:
            ERROR(funcName, "unsupported color mode");
    }

    cv::imwrite(fileName.c_str(), cvDataColor);
}

template <typename T>
void cOpenCV<T>::read(const std::string& fileName, cData3<T>& object) const
{
    const std::string   funcName("void cOpenCV::read("
                                    "const std::string& fileName, "
                                    "cData3<T>& object) const");

    cv::Mat     cvData = cv::imread(fileName.c_str(), cv::IMREAD_GRAYSCALE);

    require(cvData.data, funcName + ": could not open OpenCV file " + fileName);

    // allocate new memory
    object.memReAlloc(cVector3<size_t>((size_t) cvData.rows,
                                       (size_t) cvData.cols,
                                       1));

    #pragma omp parallel for
    for (size_t i = 0; i < object.getNelement(); i++) {
        object[i] = (T) cvData.data[i];
    }
}

template <typename T>
void cOpenCV<T>::write(const cData3<T>& object, const std::string& fileName, eNorm8bit norm8bit) const
{
    const std::string   funcName("void cOpenCV::write("
                                    "const cData3<T>& object, "
                                    "const std::string& fileName, "
                                    "eNorm8bit norm8bit) const");

    require(object.getAddrData()  != NULL, funcName + ": writing an empty object to an OpenCV file");
    require(object.getDimension() == 2,    funcName + ": writing a non-image object to an OpenCV file");

    cv::Mat     cvData((int) object.getNrow(), (int) object.getNcol(), CV_8U);
    cData3<T>   dataSave;

    switch (norm8bit) {
        case NORM8BIT_TRUE:
            dataSave = object;
            dataSave.normalize8bit();

            #pragma omp parallel for
            for(size_t i = 0; i < object.getNelement(); i++) {
                 cvData.data[i] = (uint8_t) dataSave[i];
            }
            break;
        case NORM8BIT_FALSE:
            #pragma omp parallel for
            for(size_t i = 0; i < object.getNelement(); i++) {
                 cvData.data[i] = (uint8_t) object[i];
            }
            break;
        default:
            ERROR(funcName, "unsupported norm8bit mode");
    }

    cv::imwrite(fileName.c_str(), cvData);
}

template <typename T>
void cOpenCV<T>::view2D(const cData3<T>& object, const std::string& title, eNorm8bit norm8bit)
{
    const std::string   funcName("void cOpenCV::view2D("
                                    "const cData3<T>& object, "
                                    "const std::string& title, "
                                    "eNorm8bit norm8bit)");

    require(object.getNelement()   > 0, funcName + ": viewing an empty object");
    require(object.getDimension() == 2, funcName + ": viewing a non-image object");
    require(object.getNsec()      == 1, funcName + ": viewing a non-image object");

    T         dataMin, dataMax, dataRange;
    cv::Mat   image((int) object.getNrow(), (int) object.getNcol(), CV_8U);

    switch (norm8bit) {
        case NORM8BIT_TRUE:
            dataMin   = object.getMin();
            dataMax   = object.getMax();
            dataRange = dataMax - dataMin;
            for (size_t i = 0; i < object.getNelement(); i++) {
                image.data[i] = (uint8_t) (255 * (object[i]-dataMin) / dataRange);
            }
            break;
        case NORM8BIT_FALSE:
            for (size_t i = 0; i < object.getNelement(); i++) {
                image.data[i] = (uint8_t) object[i];
            }
            break;
        default:
            ERROR("view2D", "unsupported norm8bit mode");
    }

    cv::namedWindow(title.c_str(), cv::WINDOW_NORMAL);
    cv::imshow(title.c_str(), image);

    cv::waitKey(0);
}

template <typename T>
void cOpenCV<T>::view3D(const cData3<T>& object, const std::string& title)
{
    const std::string   funcName("void cOpenCV::view3D("
                                    "const cData3<T>& object, "
                                    "const std::string& title");

    require(object.getNelement() > 0, funcName + ": viewing an empty object");
    require(object.getNcol() > 1 && object.getNsec() > 1, funcName + ": object is not 2D stack");

    size_t              nrow = 1, ncol = 1, irow, icol;
    cData3<T>           imArgNorm, imDisplay;
    unsigned int        sep = 10;

    ncol = (size_t) std::ceil(std::sqrt((double) object.getNrow()));
    nrow = (size_t) std::ceil((double) object.getNrow() / (double) ncol);

    // allocation
    imArgNorm.memAlloc(cVector3<size_t>(object.getNcol(),
                                        object.getNsec(),
                                        1));
    imDisplay.memAllocZero(cVector3<size_t>(nrow*imArgNorm.getNrow()+(nrow-1)*sep,
                                            ncol*imArgNorm.getNcol()+(ncol-1)*sep,
                                            1));

    for (size_t i = 0; i < object.getNrow(); i++) {

        imArgNorm.getSlideRow(object, i);
        imArgNorm.normalize8bit();

        irow = i / ncol;
        icol = i % ncol;

        imDisplay.opReplace(imArgNorm, cVector3<size_t>(irow*imArgNorm.getNrow()+irow*sep,
                                                        icol*imArgNorm.getNcol()+icol*sep,
                                                        0));
    }

    // display
    view2D(imDisplay, title, NORM8BIT_FALSE);
}

template <typename T>
void cOpenCV<T>::view2D(const std::string& title, unsigned int nImg, unsigned int sep, const cData3<T>* object, ...)
{
    const std::string   funcName("void cOpenCV::view2D("
                                    "const std::string& title, "
                                    "int nImg, int sep, "
                                    "const cData3<T>* object, ...)");

    require(object->getNelement()   > 0, funcName + ": viewing an empty object");
    require(object->getDimension() == 2, funcName + ": viewing a non-image object");
    require(object->getNsec()      == 1, funcName + ": viewing a non-image object");
    require(nImg > 0 && nImg <= 9,       funcName + ": improper number of objects (<= 9)");

    va_list             args;
    size_t              nrow = 1, ncol = 1, irow, icol;
    cData3<T>           *imArg, imArgNorm, imDisplay;

    switch (nImg) {
        case 1:     nrow = 1; ncol = 1; break;
        case 2:     nrow = 1; ncol = 2; break;
        case 3:     nrow = 1; ncol = 3; break;
        case 4:     nrow = 2; ncol = 2; break;
        case 5:     nrow = 2; ncol = 3; break;
        case 6:     nrow = 2; ncol = 3; break;
        case 7:     nrow = 2; ncol = 4; break;
        case 8:     nrow = 2; ncol = 4; break;
        case 9:     nrow = 3; ncol = 3; break;
        default:    ERROR("view2D", "unsupported number of objects");
    }

    // allocation
    imArgNorm.memAlloc(cVector3<size_t>(object->getNrow(),
                                        object->getNcol(),
                                        1));
    imDisplay.memAllocZero(cVector3<size_t>(nrow*object->getNrow()+(nrow-1)*sep,
                                            ncol*object->getNcol()+(ncol-1)*sep,
                                            1));

    // get the images passed
    va_start(args, object);

    // first image
    imArgNorm = *object;
    imArgNorm.normalize8bit();
    imDisplay.opReplace(imArgNorm, cVector3<size_t>(0,0,0));

    for (size_t i = 1; i < (size_t) nImg; i++) {
        imArg = va_arg(args, cData3<T>*);
        imArgNorm = *imArg;
        imArgNorm.normalize8bit();

        require(imArg->getNrow() == object->getNrow() &&
                imArg->getNcol() == object->getNcol() &&
                imArg->getNsec() == object->getNsec(),
                funcName + ": objects do not have the same size");

        irow = i / ncol;
        icol = i % ncol;

        imDisplay.opReplace(imArgNorm, cVector3<size_t>(irow*object->getNrow()+irow*sep,
                                                        icol*object->getNcol()+icol*sep,
                                                        0));
    }

    // end the number of arguments
    va_end(args);

    // display
    view2D(imDisplay, title, NORM8BIT_FALSE);
}

// instantiation
template class cOpenCV<float >;
template class cOpenCV<double>;

} // namespace gem

#endif  // __GEM_USE_OPENCV__

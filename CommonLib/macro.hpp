/***********************************************************************
 *  File:       macro.hpp
 *
 *  Purpose:    Header file for macros
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_MACRO_HPP__
#define __GEM_MACRO_HPP__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace gem {

#define GEM_PRECISON            1e-3
#define GEM_FLOATING_TOLERANCE  1e-6
#define GEM_DISTANCE_MAX        1e20

/***********************************************************************
 * enumerations
 **********************************************************************/

enum eColorName
{
    COLOR_WHITE   = 0,
    COLOR_BLACK   = 1,
    COLOR_RED     = 2,
    COLOR_GREEN   = 3,
    COLOR_BLUE    = 4,
    COLOR_CYAN    = 5,
    COLOR_MAGNETA = 6,
    COLOR_YELLOW  = 7
};

enum eNorm8bit
{
    NORM8BIT_FALSE = 0,
    NORM8BIT_TRUE  = 1
};

/***********************************************************************
 * inline functions
 **********************************************************************/

template <typename T> inline T     max(T x, T y) { return (x > y) ? x : y; }
template <typename T> inline T     min(T x, T y) { return (x < y) ? x : y; }

template <typename T> inline T     round(T x)         { return (x < 0) ? std::ceil(x-(T)0.5) : std::floor(x+(T)0.5); }
template <typename T> inline T     halfUpper(T value) { return (value % 2) ? ((value+1) / 2) : (value / 2);          }
template <typename T> inline T     halfLower(T value) { return (value % 2) ? ((value-1) / 2) : (value / 2);          }

template <typename T> inline bool  isNaN      (T value)  { return value != value; }
template <typename T> inline bool  compGreater(T x, T y) { return (x > y);        }
template <typename T> inline bool  compLesser (T x, T y) { return (x < y);        }

template <typename T> inline T     deg2rad(T deg) { return (T) deg * (T) M_PI / (T) 180 ; }
template <typename T> inline T     rad2deg(T rad) { return (T) rad * (T) 180  / (T) M_PI; }

template <typename T> inline T     pow2(T x) { return x*x;         }
template <typename T> inline T     pow3(T x) { return x*x*x;       }
template <typename T> inline T     pow4(T x) { return x*x*x*x;     }
template <typename T> inline T     pow5(T x) { return x*x*x*x*x;   }
template <typename T> inline T     pow6(T x) { return x*x*x*x*x*x; }

template <typename T> inline T     distSqr(T x1, T x2)                         { return pow2(x1-x2);                             }
template <typename T> inline T     distSqr(T x1, T y1, T x2, T y2)             { return pow2(x1-x2) + pow2(y1-y2);               }
template <typename T> inline T     distSqr(T x1, T y1, T z1, T x2, T y2, T z2) { return pow2(x1-x2) + pow2(y1-y2) + pow2(z1-z2); }

template <typename T> inline T     dist(T x1, T x2)                         { return (T) std::sqrt(distSqr(x1,x2));             }
template <typename T> inline T     dist(T x1, T y1, T x2, T y2)             { return (T) std::sqrt(distSqr(x1,y1,x2,y2));       }
template <typename T> inline T     dist(T x1, T y1, T z1, T x2, T y2, T z2) { return (T) std::sqrt(distSqr(x1,y1,z1,x2,y2,z2)); }

template <typename T> inline T     factorial(T n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n; }

inline bool                     isPow2(size_t x)                 { return ((x & (x-1)) == 0);                                    }
inline unsigned int             iDivUp(size_t a, unsigned int b) { return (unsigned int) ((a % b != 0) ? (a / b + 1) : (a / b)); }

inline void                     mpiMasterPrint(int mpiID, std::string message) { if (mpiID == 0) std::cout << message; }

template <typename T> inline
std::string num2str(T value, unsigned int precision = 0)
{
    std::ostringstream    converter;
    converter << std::fixed << std::setprecision(precision) << value;
    return converter.str();
}

/***********************************************************************
 * indexing
 **********************************************************************/

inline size_t   sub2ind(size_t iRow, size_t iCol, size_t nCol)                           { return iRow*nCol+iCol;                }
inline size_t   sub2ind(size_t iRow, size_t iCol, size_t iSec, size_t nCol, size_t nSec) { return iRow*nCol*nSec+iCol*nSec+iSec; }

template <typename T>
inline void ind2sub(size_t idx, size_t nCol, T& iRow, T& iCol)
{
    iRow = (T) std::floor((T) idx / (T) nCol);
    iCol = (T) idx - iRow * (T) nCol;
}

template <typename T>
inline void ind2sub(size_t idx, size_t nCol, size_t nSec, T& iRow, T& iCol, T& iSec)
{
    iRow = (T) std::floor((T) idx / (T) (nCol*nSec));
    iCol = (T) std::floor(((T) idx - iRow * (T) (nCol*nSec)) / (T) nSec);
    iSec = (T) idx - iRow * (T) (nCol*nSec) - iCol * (T) nSec;
}

/***********************************************************************
 * condition checking
 **********************************************************************/

inline void assure(std::ifstream &in, const std::string& filename = "")
{
    if (!(in.is_open() && in.good())) {
        fprintf(stderr, "Could not open file %s\n", filename.c_str());
        exit(EXIT_FAILURE);
    }
}

inline void assure(std::ofstream &out, const std::string& filename = "")
{
    if (!(out.is_open() && out.good())) {
        fprintf(stderr, "Could not open file %s\n", filename.c_str());
        exit(EXIT_FAILURE);
    }
}

inline void require(bool requirement, const std::string& msg = "")
{
    if (!requirement) {
        fputs("\n", stderr);
        fputs("Requirement failed: ", stderr);
        fputs(msg.c_str(), stderr);
        fputs("\n\n", stderr);
        exit(EXIT_FAILURE);
    }
}

inline void warning(bool requirement, const std::string& msg = "")
{
    if (!requirement) {
        fputs("\n", stdout);
        fputs("Warning: ", stdout);
        fputs(msg.c_str(), stdout);
        fputs("\n\n", stdout);
    }
}

/***********************************************************************
 * macros
 **********************************************************************/

#define PRINT(x)    std::cout << #x << " = " << (x) << std::endl;
#define MESSAGE(x)  std::cout << x << std::endl;

#define ERROR(function, message) \
{ \
    std::cout << std::endl << "Fatal error in " << __FILE__ \
              << " at line " << __LINE__ << "\n" \
              << function << ": " \
              << message \
              << std::endl << std::endl; \
    exit(EXIT_FAILURE); \
}

#define WARNING(function, message) \
{ \
    std::cout << std::endl << "Warning in " << __FILE__ \
              << " at line " << __LINE__ << "\n" \
              << function << ": " \
              << message \
              << std::endl << std::endl; \
}

#define CUDA_SAFE_CALL(err) \
{ \
    if (err != cudaSuccess) { \
        printf("CUDA runtime API error: %s in %s at line %d\n", \
               cudaGetErrorString(err), \
               __FILE__, __LINE__ ); \
        exit(EXIT_FAILURE); \
    } \
}

#define CUFFT_SAFE_CALL(err) \
{ \
    if (err != CUFFT_SUCCESS) { \
        printf("CUFFT error: in %s at line %d\n", \
               __FILE__, __LINE__ ); \
        exit(EXIT_FAILURE); \
    } \
}

#define CUDA_ARRAY_PRINT(array, length, T) \
{ \
    T   *arrayTmp = NULL; \
    array_new(arrayTmp, length); \
    cuda_array_memcpy_d2h(arrayTmp, array, length); \
    array_print(arrayTmp, length); \
    array_delete(arrayTmp); \
}

} // namespace gem

#endif

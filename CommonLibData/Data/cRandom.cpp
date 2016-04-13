/***********************************************************************
 *  File:       cRandom.cpp
 *
 *  Purpose:    Implementation of a random class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cRandom.hpp"

#ifdef __GEM_USE_STD11__
#include <random>
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#pragma GCC diagnostic pop
#endif

namespace gem {

template <typename T>
void cRandom<T>::randn(cData3<T>& data, T mean, T std, const cSize3& size)
{
    const std::string   funcName("void cRandom<T>::randn(cData3<T>& data, "
                                    "T mean, T std, const cSize3& size)");

    require(size.getProduct() > 0, funcName + ": an empty object is requested");

    data.memReAlloc(size);

#ifdef __GEM_USE_STD11__
    std::random_device              ranDev;
    std::mt19937_64                 ranGen(ranDev());
    std::normal_distribution<T>     ranDist(mean, std);
#else
    boost::random::random_device            ranDev;
    boost::mt19937_64                       ranGen(ranDev());
    boost::random::normal_distribution<T>   ranDist(mean, std);
#endif

    for (size_t i = 0; i < data.getNelement(); i++) {
        data[i] = ranDist(ranGen);
    }
}

template <typename T>
void cRandom<T>::rand(cData3<T>& data, T min, T max, const cSize3& size)
{
    const std::string   funcName("void cRandom<T>::rand(cData3<T>& data, "
                                    "T min, T max, const cSize3& size)");

    require(min <= max,            funcName + ": min > max");
    require(size.getProduct() > 0, funcName + ": an empty object is requested");

    data.memReAlloc(size);

#ifdef __GEM_USE_STD11__
    std::random_device                  ranDev;
    std::mt19937_64                     ranGen(ranDev());
    std::uniform_real_distribution<T>   ranDist(min, max);
#else
    boost::random::random_device                    ranDev;
    boost::mt19937_64                               ranGen(ranDev());
    boost::random::uniform_real_distribution<T>     ranDist(min, max);
#endif

    for (size_t i = 0; i < data.getNelement(); i++) {
        data[i] = ranDist(ranGen);
    }
}

template <typename T>
void cRandom<T>::randZero(cData3<T>& data, T min, T max, T threshold)
{
    const std::string   funcName("void cRandom<T>::randZero(cData3<T>& data, "
                                    "T min, T max, T threshold)");

    data.requireNonEmpty(funcName);

    cData3<T>       noise;

    rand(noise, min, max, data.getSize());

    array_remove_rezo(noise.getAddrData(), data.getNelement(),
                       data.getAddrData(),
                      threshold);
}

#ifdef __GEM_USE_CUDA__
template <typename T>
void cRandom<T>::randn(cuData3<T>& data, T mean, T std, const cSize3& size, size_t seed)
{
    const std::string   funcName("void cRandom<T>::randn(cuData3<T>& data, "
                                    "T mean, T std, const cSize3& size, size_t seed)");

    require(size.getProduct() > 0, funcName + ": an empty object is requested");

    data.memReAlloc(size);

    cuda_array_random_normal(data.getAddrData(),
                             data.getNelement(),
                             mean, std, seed);
}
#endif

#ifdef __GEM_USE_CUDA__
template <typename T>
void cRandom<T>::rand(cuData3<T>& data, T min, T max, const cSize3& size, size_t seed)
{
    const std::string   funcName("void cRandom<T>::rand(cuData3<T>& data, "
                                    "T min, T max, const cSize3& size, size_t seed)");

    require(min <= max,            funcName + ": min > max");
    require(size.getProduct() > 0, funcName + ": an empty object is requested");

    data.memReAlloc(size);

    cuda_array_random_uniform(data.getAddrData(),
                              data.getNelement(),
                              min, max, seed);
}
#endif

#ifdef __GEM_USE_CUDA__
template <typename T>
void cRandom<T>::randZero(cuData3<T>& data, T min, T max, T threshold)
{
    const std::string   funcName("void cRandom<T>::randZero(cuData3<T>& data, "
                                    "T min, T max, T threshold)");

    data.requireNonEmpty(funcName);

    cuData3<T>      noise;

    rand(noise, min, max, data.getSize());

    cuda_array_remove_rezo(noise.getAddrData(), data.getNelement(),
                            data.getAddrData(),
                           threshold);
}
#endif

// instantiation
template class cRandom<float >;
template class cRandom<double>;

} // namespace gem

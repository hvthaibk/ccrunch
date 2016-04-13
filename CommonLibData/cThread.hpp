/***********************************************************************
 *  File:       cThread.hpp
 *
 *  Purpose:    Header file for thread-safe data structures
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CTHREAD_HPP__
#define __GEM_CTHREAD_HPP__

#include "macro.hpp"

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>

namespace gem {

template<class T>
class cQueueConcurrent
{
private:
    std::queue<std::shared_ptr<T>>    queueHandle;
    mutable std::mutex                mutexHandle;
    std::condition_variable           condFlag;
    size_t                            sizeMax;

public:
     cQueueConcurrent(unsigned int length = 100) : sizeMax(length) {}
    ~cQueueConcurrent() { require(empty(), "cQueueConcurrent is not empty at destruction"); }

    // for message queue
    void push(T const data)
    {
        std::shared_ptr<T>              element(std::make_shared<T>(std::move(data)));
        std::lock_guard<std::mutex>     lockHandle(mutexHandle);

        queueHandle.push(element);
        condFlag.notify_all();
    }

    bool try_pop(T& data)
    {
        std::lock_guard<std::mutex>     lockHandle(mutexHandle);

        if (!queueHandle.empty()) {
            data = std::move(*queueHandle.front());
            queueHandle.pop();
            condFlag.notify_one();

            return true;
        }
        else {
            return false;
        }
    }

    // for data queue
    void wait_and_push(T const data)
    {
        std::shared_ptr<T>               element(std::make_shared<T>(std::move(data)));
        std::unique_lock<std::mutex>     lockHandle(mutexHandle);

        condFlag.wait(lockHandle,
                      [this]{return (queueHandle.size() < sizeMax);});

        queueHandle.push(element);
        condFlag.notify_all();
    }

    /*bool wait_and_push(T const data, std::chrono::microseconds waitingTime)
    {
        std::shared_ptr<T>               element(std::make_shared<T>(std::move(data)));
        std::unique_lock<std::mutex>     lockHandle(mutexHandle);

        if (condFlag.wait_for(lockHandle, waitingTime,
                     [this]{return (queueHandle.size() < sizeMax);})) {
            queueHandle.push(element);
            condFlag.notify_all();
            return true;
        }
        else {
            return false;
        }
    }*/

    /*void wait_and_pop(T& data)
    {
        std::unique_lock<std::mutex>     lockHandle(mutexHandle);

        condFlag.wait(lockHandle, [this]{return !queueHandle.empty();});

        data = std::move(*queueHandle.front());
        queueHandle.pop();
        condFlag.notify_all();
    }*/

    bool wait_and_pop(T& data, std::chrono::microseconds waitingTime)
    {
        std::unique_lock<std::mutex>     lockHandle(mutexHandle);

        if (condFlag.wait_for(lockHandle, waitingTime,
                              [this]{return !queueHandle.empty();})) {
            data = std::move(*queueHandle.front());
            queueHandle.pop();
            condFlag.notify_all();

            return true;
        }
        else {
            return false;
        }
    }

    // for condition checking
    bool empty() const
    {
        std::lock_guard<std::mutex>     lockHandle(mutexHandle);

        return queueHandle.empty();
    }

    /*size_t size() const
    {
        std::lock_guard<std::mutex>     lockHandle(mutexHandle);

        return queueHandle.size();
    }*/
};

} // namespace gem

#endif

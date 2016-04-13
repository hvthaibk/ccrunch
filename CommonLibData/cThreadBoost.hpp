/***********************************************************************
 *  File:       cThreadBoost.hpp
 *
 *  Purpose:    Header file for thread-safe data structures
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CTHREADBOOST_HPP__
#define __GEM_CTHREADBOOST_HPP__

#include "config.hpp"
#include "macro.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include "boost/thread.hpp"
#pragma GCC diagnostic pop

#include <queue>

namespace gem {

template<class T>
class cQueueConcurrent
{
private:
    std::queue<T>                queueHandle;
    mutable boost::mutex         mutexHandle;
    boost::condition_variable    condFlag;
    size_t                       sizeMax;

public:
     cQueueConcurrent(unsigned int length = 100) : sizeMax(length) {}
    ~cQueueConcurrent() { require(empty(), "cQueueConcurrent is not empty at destruction"); }

    // for message queue
    void push(T data)
    {
        boost::mutex::scoped_lock     lockHandle(mutexHandle);

        queueHandle.push(data);
        lockHandle.unlock();
        condFlag.notify_all();
    }

    bool try_pop(T& data)
    {
        boost::mutex::scoped_lock     lockHandle(mutexHandle);

        if (!queueHandle.empty()) {
            data = queueHandle.front();
            queueHandle.pop();
            condFlag.notify_one();

            return true;
        }
        else {
            return false;
        }
    }

    // for data queue
    void wait_and_push(T data)
    {
        boost::mutex::scoped_lock     lockHandle(mutexHandle);

        condFlag.wait(lockHandle, isFull(&queueHandle, sizeMax));
        queueHandle.push(data);
        lockHandle.unlock();
        condFlag.notify_all();
    }

    bool wait_and_pop(T& data, boost::posix_time::microseconds waitTime)
    {
        boost::system_time            timeout = boost::get_system_time() + waitTime;
        boost::mutex::scoped_lock     lockHandle(mutexHandle);

        if (condFlag.timed_wait(lockHandle, timeout,
                                isEmpty(&queueHandle))) {
            data = queueHandle.front();
            queueHandle.pop();
            condFlag.notify_all();

            return true;
        }
        else {
            return false;
        }
    }

    bool empty() const
    {
        boost::mutex::scoped_lock     lockHandle(mutexHandle);

        return queueHandle.empty();
    }

    struct isEmpty
    {
        std::queue<T>    *queue;

        isEmpty(std::queue<T> *queue_): queue(queue_) {}

        bool operator()() const { return !queue->empty(); }
    };

    struct isFull
    {
        std::queue<T>    *queue;
        size_t           sizeMax;

        isFull(std::queue<T> *queue_,
               size_t sizeMax_) : queue(queue_),
               sizeMax(sizeMax_) {}

        bool operator()() const { return queue->size() < sizeMax; }
    };
};

} // namespace gem

#endif

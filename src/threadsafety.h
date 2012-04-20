//
//  thread_safety.h
//  express
//
//  Created by Adam Roberts on 1/25/12.
//  Copyright 2012 Adam Roberts. All rights reserved.
//

#ifndef express_thread_safety_h
#define express_thread_safety_h

#include <boost/thread.hpp>
#include <queue.h>

class Fragment;

class ThreadSafeFragQueue
{
    std::queue<Fragment*> _queue;
    size_t _max_size;
    boost::mutex _mut;
    boost::condition_variable _cond;
    
public:
    ThreadSafeFragQueue(size_t max_size);
    Fragment* pop(bool block=true);
    void push(Fragment* frag);
    bool empty(bool block=false);
};

/**
 * The ParseThreadSafety struct stores objects to allow for parsing to safely occur
 * on a separate thread from processing.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
struct ParseThreadSafety
{
    ThreadSafeFragQueue proc_in;
    ThreadSafeFragQueue proc_on;
    ThreadSafeFragQueue proc_out;
    ParseThreadSafety(size_t q_size) : proc_in(q_size), proc_on(q_size), proc_out(q_size) {}
};

#endif

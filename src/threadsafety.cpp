//
//  threadsafety.cpp
//  express
//
//  Created by Adam Roberts on 4/19/12.
//  Copyright (c) 2012 Adam Roberts. All rights reserved.
//

#include "threadsafety.h"
#include "fragments.h"

ThreadSafeFragQueue::ThreadSafeFragQueue(size_t max_size)
    : _max_size(max_size) {
}

Fragment* ThreadSafeFragQueue::pop(bool block) {
  boost::unique_lock<boost::mutex> lock(_mut);
  while (_queue.empty()) {
    if (!block) {
            return NULL;
    }
        _cond.wait(lock);
  }
    
  _cond.notify_all();
  Fragment* res = _queue.front();
  _queue.pop();
  return res;
}

void ThreadSafeFragQueue::push(Fragment* frag) {
  boost::unique_lock<boost::mutex> lock(_mut);
  while (_queue.size() == _max_size) {
    _cond.wait(lock);
  }
    
  _cond.notify_all();
  return _queue.push(frag);
}

bool ThreadSafeFragQueue::is_empty(bool block) {
  boost::unique_lock<boost::mutex> lock(_mut);
  while (!_queue.empty()) {
    if (!block) {
      return false;
    }
    _cond.wait(lock);
  }
  return true;
}

//
//  robertsfilter.cpp
//  express
//
//  Created by Adam Roberts on 5/1/12.
//  Copyright (c) 2012 UC Berkeley. All rights reserved.
//

#include "robertsfilter.h"

using namespace std;

RobertsFilter::RobertsFilter(size_t local_size, size_t global_size)
    : _global_vector(global_size),
      _local_size(local_size),
      _global_size(global_size) {
}

bool RobertsFilter::test_and_push(const string& key) {
  if (_local_set.count(key) || _global_set.count(key)) {
    return true;
  }
    
  _local_set.insert(key);
  _local_queue.push(key);
  if (_local_set.size() > _local_size) {
    size_t r = _global_set.size();
    if (_global_set.size() == _global_size) {
      r = (size_t)(rand()/double(RAND_MAX)*(_global_size-1));
      _global_set.erase(_global_vector[r]);
    }

    _local_set.erase(_local_queue.front());
    _global_set.insert(_local_queue.front());
    _global_vector[r] = _local_queue.front();
    _local_queue.pop();
  }
  return false;
}

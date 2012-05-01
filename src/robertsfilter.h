//
//  robertsfilter.h
//  express
//
//  Created by Adam Roberts on 5/1/12.
//  Copyright (c) 2012 Adam Roberts. All rights reserved.
//

#ifndef express_robertsfilter_h
#define express_robertsfilter_h

#include <queue>
#include <vector>
#include <boost/unordered_set.hpp>

static size_t DEFAULT_LOC_SIZE = 10000;
static size_t DEFAULT_GLOB_SIZE = 100000;

class RobertsFilter
{
    std::queue<std::string> _local_queue;
    boost::unordered_set<std::string> _local_set;
    std::vector<std::string> _global_vector;
    boost::unordered_set<std::string> _global_set;
    
    size_t _local_size;
    size_t _global_size;

public:
    RobertsFilter(size_t local_size=DEFAULT_LOC_SIZE, size_t global_size=DEFAULT_GLOB_SIZE);
    bool test_and_push(const std::string& frag_name);
};


#endif

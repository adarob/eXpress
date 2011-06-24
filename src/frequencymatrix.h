#ifndef FREQUENCYMATRIX_H
#define FREQUENCYMATRIX_H

//
//  frequencymatrix.h
//  expressionline2
//
//  Created by Adam Roberts on 4/23/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
#include "boost/multi_array.hpp"
#include <vector>

class FrequencyMatrix
{
    std::vector<double> _array;
    std::vector<double> _colsums;
    size_t _M;
    size_t _N;
    
public:
    FrequencyMatrix(){};
    FrequencyMatrix(size_t m, size_t n, double alpha);
    
    double operator()(size_t i, size_t j) const;
    double operator()(size_t i) const;
    void increment(size_t i, size_t j, double incr_amt);
    void increment(size_t i, double incr_amt);
};

#endif
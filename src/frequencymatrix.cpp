//
//  frequencymatrix.cpp
//  express
//
//  Created by Adam Roberts on 4/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include <assert.h>
#include "frequencymatrix.h"
#include "main.h"
#include <boost/array.hpp>

using namespace std;

FrequencyMatrix::FrequencyMatrix(size_t m, size_t n, double alpha)
: _array(m*n, log(alpha)),
  _rowsums(m,log(n*alpha)),
  _M(m),
  _N(n)
{}


double FrequencyMatrix::operator()(size_t i, size_t j) const
{
    assert(i*_N+j < _M*_N);
    assert(!isnan(_array[i*_N+j]-_rowsums[i]));
    return _array[i*_N+j]-_rowsums[i];
}

double FrequencyMatrix::operator()(size_t k) const
{
    return operator()(0, k);
}


void FrequencyMatrix::increment(size_t i, size_t j, double incr_amt)
{
    assert(i*_N+j < _M*_N);
    _array[i*_N+j] = log_sum(_array[i*_N+j], incr_amt);
    _rowsums[i] = log_sum(_rowsums[i], incr_amt);
    assert(!isnan(_rowsums[i]) && !isinf(_rowsums[i]));
}

void FrequencyMatrix::increment(size_t k, double incr_amt)
{
    increment(0, k, incr_amt);
}

//
//  frequencymatrix.cpp
//  expressionline2
//
//  Created by Adam Roberts on 4/23/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <assert.h>
#include "frequencymatrix.h"
#include <boost/array.hpp>

using namespace std;

FrequencyMatrix::FrequencyMatrix(size_t m, size_t n, double alpha)
: _array(m*n, alpha),
  _colsums(m,n*alpha),
  _M(m),
  _N(n)
{}


double FrequencyMatrix::operator()(size_t i, size_t j) const
{
    assert(i*_N+j < _M*_N);
    assert(_array[i*_N+j] > 0);
    assert(_colsums[i] > 0);
    assert(!isnan(_array[i*_N+j]) && !isinf(_array[i*_N+j]));
    return _array[i*_N+j]/_colsums[i];
}

double FrequencyMatrix::operator()(size_t i) const
{
    return operator()(0, i);
}

void FrequencyMatrix::increment(size_t i, size_t j, double incr_amt)
{
    assert(i*_N+j < _M*_N);
    assert(!isnan(incr_amt) && !isinf(incr_amt));
    _array[i*_N+j] += incr_amt;
    _colsums[i] += incr_amt;
}

void FrequencyMatrix::increment(size_t i, double incr_amt)
{
    increment(0, i, incr_amt);
}

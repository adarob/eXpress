//
//  frequencymatrix.cpp
//  express
//
//  Created by Adam Roberts on 4/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include <cassert>
#include "frequencymatrix.h"
#include "main.h"
#include <boost/array.hpp>

using namespace std;

FrequencyMatrix::FrequencyMatrix(size_t m, size_t n, double alpha, bool logged)
: _array(m*n, logged ? log(alpha):alpha),
  _rowsums(m, logged ? log(n*alpha):n*alpha),
  _M(m),
  _N(n),
  _logged(logged)
{}


double FrequencyMatrix::operator()(size_t i, size_t j) const
{
    assert(i*_N+j < _M*_N);
    if (_logged)
        return _array[i*_N+j]-_rowsums[i];
    else
        return _array[i*_N+j]/_rowsums[i];
}

double FrequencyMatrix::operator()(size_t k) const
{
    return operator()(0, k);
}


void FrequencyMatrix::increment(size_t i, size_t j, double incr_amt)
{
    size_t k = i*_N+j;
    assert(k < _M*_N);
    if (_logged)
    {
        _array[k] = log_sum(_array[k], incr_amt);
        _rowsums[i] = log_sum(_rowsums[i], incr_amt);
    }
    else
    {
        _array[k] += incr_amt;
        _rowsums[i] += incr_amt;
    }
//    assert(!isnan(_rowsums[i]) && !isinf(_rowsums[i]));
}

void FrequencyMatrix::increment(size_t k, double incr_amt)
{
    increment(0, k, incr_amt);
}

void FrequencyMatrix::set_logged(bool logged)
{
    if (logged == _logged)
        return;
    
    if (logged)
    {
        for(size_t i = 0; i < _M*_N; ++i)
        {
            _array[i] = log(_array[i]);
        }
        for(size_t i = 0; i < _M; ++i)
        {
            _rowsums[i] = log(_rowsums[i]);
        }
    }
    else
    {
        for(size_t i = 0; i < _M*_N; ++i)
        {
            _array[i] = sexp(_array[i]);
        }
        for(size_t i = 0; i < _M; ++i)
        {
            _rowsums[i] = sexp(_rowsums[i]);
        }
    }
    
    _logged = logged;
    
}
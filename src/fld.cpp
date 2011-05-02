//
//  fld.cpp
//  expressionline
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "fld.h"
#include <numeric>

using namespace std;

FLD::FLD(double alpha, size_t max_val) : _min(max_val)
{
    _hist = vector<double>(max_val+1, alpha);
    _total = max_val * alpha;
}

size_t FLD::max_val() const
{
    return _hist.size()-1;
}

void FLD::add_val(size_t len, double mass)
{
    if (len > max_val()) return;
    if (len < _min) _min = len;
    _hist[len] += mass;
    _total += mass;
}

double FLD::pdf(size_t len) const
{
    if (len > max_val())
        return 0.0;
    return _hist[len]/_total;
}

double FLD::cdf(size_t len) const 
{
    return accumulate(_hist.begin(),_hist.begin()+len+1,0.0)/_total;
}

double FLD::npdf(size_t len, size_t max) const
{
    return _hist[len]/accumulate(_hist.begin(),_hist.begin()+max+1,0.0);
}

double FLD::total() const
{
    return _total;
}
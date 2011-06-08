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
    _num_obs = max_val * alpha;
    _sum = (max_val)*(max_val+1)*alpha/2;
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
    _num_obs += mass;
    _sum += mass*len;
}

double FLD::pdf(size_t len) const
{
    if (len > max_val())
        return 0.0;
    return _hist[len]/_num_obs;
}

double FLD::cdf(size_t len) const 
{
    return accumulate(_hist.begin(),_hist.begin()+len+1,0.0)/_num_obs;
}

double FLD::npdf(size_t len, size_t max) const
{
    return _hist[len]/accumulate(_hist.begin(),_hist.begin()+max+1,0.0);
}

double FLD::num_obs() const
{
    return _num_obs;
}

double FLD::mean() const
{
    return _sum/num_obs();
}
//
//  fld.cpp
//  expressionline
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "fld.h"
#include <numeric>
#include <boost/assign.hpp>
#include <iostream>
#include <fstream>

using namespace std;

const vector<double> KERNEL = boost::assign::list_of(0.0625)(0.25)(0.375)(0.25)(0.0625); 

FLD::FLD(double alpha, size_t max_val) : 
    _hist(max_val+1, alpha),
    _num_obs(max_val * alpha), 
    _sum((max_val)*(max_val+1)*alpha/2),
    _min(max_val)
{
    assert(KERNEL.size() % 2 == 1);
}

size_t FLD::max_val() const
{
    return _hist.size()-1;
}

void FLD::add_val(size_t len, double mass)
{
    if (len > max_val()) return;
    if (len < _min) _min = len;
    
    size_t offset = len - KERNEL.size()/2; 
    
    for (size_t i = 0; i < KERNEL.size(); i++)
    {
        double k_mass = mass * KERNEL[i];
        _hist[offset] += k_mass;
        _sum += (offset++)*k_mass;
    }
    _num_obs += mass;
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

void FLD::output(string path) const
{
    string filename = path + "/frag_len_hist.tab";
    ofstream outfile(filename.c_str());
    outfile<<"Length\tMass\n";
    for(size_t i = 0; i < max_val()+1; i++)
    {
        outfile << i << '\t' << _hist[i]; 
    }
    outfile.close();
}
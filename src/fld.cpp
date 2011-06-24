//
//  fld.cpp
//  expressionline
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "fld.h"
#include "main.h"
#include <numeric>
#include <boost/assign.hpp>
#include <iostream>
#include <fstream>

using namespace std;

const vector<double> KERNEL = boost::assign::list_of(-1.20411998)(-0.602059991)(-0.425968732)(-0.602059991)(-1.20411998); 

FLD::FLD(double alpha, size_t max_val) : 
    _hist(max_val+1, log(alpha)),
    _num_obs(log(max_val * alpha)), 
    _sum(log((max_val)*(max_val+1)*alpha/2)),
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
        double k_mass = mass + KERNEL[i];
        _hist[offset] = log_sum(_hist[offset], k_mass);
        _sum = log_sum(_sum, log(offset++)+k_mass);
    }
    _num_obs = log_sum(_num_obs, mass);
}

double FLD::pdf(size_t len) const
{
    if (len > max_val())
        return HUGE_VAL;
    return _hist[len]-_num_obs;
}

double FLD::num_obs() const
{
    return _num_obs;
}

double FLD::mean() const
{
    return _sum - num_obs();
}

void FLD::output(string path) const
{
    string filename = path + "/frag_len_hist.tab";
    ofstream outfile(filename.c_str());
    outfile<<"Length\tMass\n";
    for(size_t i = 0; i < max_val()+1; i++)
    {
      outfile << i << '\t' << sexp(_hist[i]) << '\n'; 
    }
    outfile.close();
}

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
    _tot_mass(log(max_val * alpha)), 
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
    _tot_mass = log_sum(_tot_mass, mass);
}

double FLD::pdf(size_t len) const
{
    if (len > max_val())
        return HUGE_VAL;
    return _hist[len]-_tot_mass;
}

double FLD::tot_mass() const
{
    return _tot_mass;
}

double FLD::mean() const
{
    return _sum - tot_mass();
}

string FLD::to_string() const
{
    string s = "";
    char buffer[50];
    for(size_t i = 0; i < max_val()+1; i++)
    {
        sprintf(buffer, "%e,",sexp(_hist[i]));
        s += buffer; 
    }
    s.erase(s.length()-1,1);
    return s;
}

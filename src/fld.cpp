//
//  fld.cpp
//  express
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
#include <boost/math/distributions/normal.hpp>

using namespace std;

const vector<double> KERNEL = boost::assign::list_of(-2.7725887222397811)(-1.3862943611198906)(-0.98082925301172619)(-1.3862943611198906)(-2.7725887222397811); 

FLD::FLD(double alpha, size_t max_val, size_t mean, size_t std_dev) : 
    _hist(max_val+1),
    _tot_mass(HUGE_VAL), 
    _sum(HUGE_VAL),
    _min(max_val)
{
    assert(KERNEL.size() % 2 == 1);
    boost::math::normal norm(mean, std_dev);
    
    double tot = log(alpha*max_val);
    
    for (size_t i = 1; i <= max_val; ++i)
    {  
        double mass = tot + log(boost::math::cdf(norm,i+0.5) - boost::math::cdf(norm,i-0.5));
        _hist[i] = mass;
        _sum = log_sum(_sum, log((double)i)+mass);
        _tot_mass = log_sum(_tot_mass, mass);
    }
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
        _sum = log_sum(_sum, log((double)offset++)+k_mass);
    }
    _tot_mass = log_sum(_tot_mass, mass);
}

double FLD::pdf(size_t len) const
{
    if (len < 1 || len > max_val())
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
    for(size_t i = 0; i <= max_val(); i++)
    {
        sprintf(buffer, "%e\t",sexp(pdf(i)));
        s += buffer; 
    }
    s.erase(s.length()-1,1);
    return s;
}

void FLD::append_output(ofstream& outfile) const
{
    outfile << ">Fragment Length Distribution (0-" << max_val() << ")\n";
    outfile << to_string() << endl;
}

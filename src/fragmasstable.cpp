//
//  fragmasstable.cpp
//  express
//
//  Created by Adam Roberts on 9/7/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "fragmasstable.h"
#include "main.h"

size_t FragMassTable::n_to_index (size_t n)
{
    if (n < 8)
        return n;
    
    size_t msb = INT_BIT - __builtin_clz((int)n);
    n = n >> (msb-3);
    return 4*(msb-2) + (n & 3);
}

FragMassTable::FragMassTable(double ff_param)
: _ff_param(ff_param),
  _n(0),
  _mass_n(0.0),
  _cum_mass_n(HUGE_VAL)
{}

double FragMassTable::next_frag_mass()
{
    if (_n > 0)
    {
        _cum_mass_n = log_sum(_cum_mass_n, _mass_n);
        _mass_n += _ff_param*log((double)_n) - log(pow(_n+1,_ff_param) - 1);
    }

    if (_mass_table.size() == n_to_index(_n)) 
    {
        _n_table.push_back(_n);
        _mass_table.push_back(_mass_n);
        _cum_mass_table.push_back(_cum_mass_n);
    }
    
    ++_n;
    
    return _mass_n;
}

double FragMassTable::next_frag_mass(size_t n, double curr_mass)
{
    if (n)
    {
        return curr_mass +_ff_param*log((double)n) - log(pow(n+1,_ff_param) - 1);
    }
    return 0.0;
}

size_t FragMassTable::nearest_stored_mass(size_t n, double& mass, double& cum_mass)
{
    size_t i = n_to_index(n);
    mass = _mass_table[i];
    cum_mass = _cum_mass_table[i];
    return _n_table[i];
}

//
//  markovmodel.cpp
//  express
//
//  Created by Adam Roberts on 1/24/12.
//  Copyright 2012 Adam Roberts. All rights reserved.
//

#include "markovmodel.h"
#include "frequencymatrix.h"
#include "sequence.h"
#include "main.h"

using namespace std;

const double LOG_QUARTER = log(0.25);

MarkovModel::MarkovModel(size_t order, size_t window_size, size_t num_pos, double alpha, bool logged)
: _order((int)order),
  _window_size((int)window_size),
  _num_pos((int)num_pos),
  _params(num_pos, FrequencyMatrix(pow((double)NUM_NUCS, (double)order), NUM_NUCS, alpha, logged)),
  _logged(logged)
{}


void MarkovModel::update(const Sequence& seq, int left, double mass)
{
    int i = 0;
    int j = left;
    int seq_len = (int)seq.length();
    
    size_t cond = 0;
    
    if (left < _order)
    {
        i = _order-left;
        for (j=0; j < min(_order, i); j++)
        {
            cond = (cond << 2) + seq[j];
        }
    }
    
    while(i < _window_size && j < seq_len)
    {
        size_t curr = seq[j];
        _params[min(i,_num_pos-1)].increment(cond, curr, mass);
        cond = (cond << 2) + curr;
        if (i >= _order)
            cond -= seq[j-_order] << (_order*2);
        i++;
        j++;
    }
}

void MarkovModel::fast_learn(const Sequence& seq, double mass, const vector<double>& fl_cdf)
{
    assert(_num_pos==_order+1);
    if (seq.length() < _order)
        return;
    
    size_t cond = 0;
    for (int i = 0; i < _order; ++i)
    {
        cond = (cond << 2) + seq[i];
    }
    
    for(size_t i = _order; i < seq.length(); ++i)
    {
        size_t curr = seq[i];
        
        double mass_i = mass;
        if (seq.length()-i < fl_cdf.size())
        {
            mass_i += fl_cdf[seq.length()-1];
        }
        
        _params[_order].increment(cond, curr, mass_i);
        cond = (cond << 2) + curr;
        cond -= seq[i-_order] << (_order*2);
    }
}

void MarkovModel::calc_marginals()
{
    assert(_num_pos==_order+1);
    for (int i = 0; i < _order; ++i)
    {
        for (size_t cond = 0; cond < pow((double)NUM_NUCS, (double)(_order)); cond++)
        {
            for(size_t nuc = 0; nuc < NUM_NUCS; ++nuc)
            {
                _params[i].increment(cond & ((1 << (2*i))-1) , nuc, _params[_order].arr((cond << 2) + nuc));
            }
        }
    }
}

double MarkovModel::transition_prob(size_t p, size_t cond, size_t curr) const
{
    assert(p < _params.size());
    return _params[p](cond, curr);
}

double MarkovModel::seq_prob(const Sequence& seq, int left) const
{

    int i = 0;
    int j = left;
    int seq_len = (int)seq.length();
    
    size_t cond = 0;
    double v = (_logged) ? 0:1;

    if (left < _order)
    {
        i = _order-left;
        for (j=0; j < min(_order, i); j++)
        {
            cond = (cond << 2) + seq[j];
        }
        if (_logged)
            v = i*LOG_QUARTER;
        else
            v = pow(0.25, (double)i);
    }
        
    while(i < _window_size && j < seq_len)
    {
        size_t curr = seq[j];
        if (_logged)
            v += _params[min(i,_num_pos-1)](cond, curr);
        else
            v *= _params[min(i,_num_pos-1)](cond, curr);
        cond = (cond << 2) + curr;
        if (i >= _order)
            cond -= seq[j-_order] << (_order*2);
        i++;
        j++;
    }
    
    return v;
}

void MarkovModel::set_logged(bool logged)
{
    if (logged == _logged)
        return;
    _logged = logged;
    for (size_t i = 0; i < _params.size(); ++i)
    {
        _params[i].set_logged(logged);
    }
}

double MarkovModel::marginal_prob(size_t w, size_t nuc) const
{
    assert(_logged);
    assert(w < _params.size());
    double marg = HUGE_VAL;
    double tot = HUGE_VAL;
    for (size_t cond = 0; cond < pow((double)NUM_NUCS, (double)(_order)); cond++)
    {
        marg = log_sum(marg, _params[w].arr((cond << 2) + nuc));
        tot = log_sum(tot, _params[w].row(cond));
    }
    return marg-tot;
}

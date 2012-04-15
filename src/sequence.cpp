//
//  sequence.cpp
//  express
//
//  Created by Adam Roberts on 1/24/12.
//  Copyright 2012 Adam Roberts. All rights reserved.
//

#include "sequence.h"
#include <cassert>
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;


SequenceFwd::SequenceFwd():  _ref_seq(NULL), _len(0) {}

SequenceFwd::SequenceFwd(const std::string& seq, bool rev, bool prob) : _ref_seq(NULL), _len(seq.length()), _prob(prob)
{
    if (prob)
    {
        _est_seq = FrequencyMatrix<float>(seq.length(), NUM_NUCS, 0.01);
        _obs_seq = FrequencyMatrix<float>(seq.length(), NUM_NUCS, HUGE_VAL);
        _exp_seq = FrequencyMatrix<float>(seq.length(), NUM_NUCS, HUGE_VAL);
    }
    set(seq, rev);
}

SequenceFwd::SequenceFwd(const SequenceFwd& other) : _ref_seq(NULL), _obs_seq(other._obs_seq), _exp_seq(other._exp_seq), _len(other.length()), _prob(other._prob)
{
    if (other._ref_seq)
    {
        char* ref_seq = new char[_len];
        std::copy(other._ref_seq, other._ref_seq + _len, ref_seq);
        _ref_seq = ref_seq;
    }
}

SequenceFwd& SequenceFwd::operator=(const SequenceFwd& other)
{
    if (other._ref_seq)
    {
        _len = other.length();
        char* ref_seq = new char[_len];
        std::copy(other._ref_seq, other._ref_seq + _len, ref_seq);
        _ref_seq = ref_seq;
        _obs_seq = other._obs_seq;
        _exp_seq = other._exp_seq;
        _prob = other._prob;
    }
    return *this;
}


SequenceFwd::~SequenceFwd()
{
    if (_ref_seq)
        delete _ref_seq;
}


void SequenceFwd::set(const std::string& seq, bool rev)
{
    if (_ref_seq)
        delete _ref_seq;
    
    char* ref_seq = new char[seq.length()]; 
    for(size_t i = 0; i < seq.length(); i++)
    {
        ref_seq[i] = (rev) ? complement(ctoi(seq[seq.length()-1-i])) : ctoi(seq[i]);
        if (_prob)
        {
            _est_seq.increment(i, ref_seq[i], log(2));
        }
    }
    _ref_seq = ref_seq;
    _len = seq.length();
}

size_t SequenceFwd::operator[](const size_t index) const
{
    assert(index < _len);
    if (_prob)
    {
        return _est_seq.mode(index);
    }
    return _ref_seq[index]; 
}

size_t SequenceFwd::get_ref(const size_t index) const
{
    assert(index < _len);
    assert(_ref_seq[index] == operator[](index));
    return _ref_seq[index]; 
}

float SequenceFwd::get_prob(const size_t index, const size_t nuc) const
{
    assert(_prob);
    return _est_seq(index, nuc);
}

void SequenceFwd::update_est(const size_t index, const size_t nuc, float mass)
{
    assert(_prob);
    _est_seq.increment(index, nuc, mass);
}

void SequenceFwd::update_obs(const size_t index, const size_t nuc, float mass)
{
    assert(_prob);
    _obs_seq.increment(index, nuc, mass);
}

void SequenceFwd::update_exp(const size_t index, const size_t nuc, float mass)
{
    assert(_prob);
    _exp_seq.increment(index, nuc, mass);
}

void SequenceFwd::calc_p_vals(vector<double>& p_vals) const
{
    p_vals = vector<double>(_len, 1.0);
    boost::math::chi_squared_distribution<double> chisq(3);
    double obs_n,exp_n;
    for (size_t i = 0; i < _len; ++i)
    {
        if (_obs_seq.total(i)==HUGE_VAL)
            continue;
        
        double S = 0;
        for (size_t nuc = 0; nuc < NUM_NUCS; ++nuc)
        {
            obs_n = sexp(_obs_seq(i,nuc,false));
            exp_n = sexp(_exp_seq(i,nuc,false));
            S += (obs_n-exp_n)*(obs_n-exp_n)/exp_n;
        }
        
        p_vals[i] -= boost::math::cdf(chisq,S);
    }
}

void SequenceRev::calc_p_vals(vector<double>& p_vals) const
{
    vector<double> temp;
    _seq->calc_p_vals(temp);
    p_vals = vector<double>(length());
    for(size_t i = 0; i < length(); i++)
    {
        p_vals[i] = temp[length()-i-1];
    }
}
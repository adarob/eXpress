/*
 *  biascorrection.h
 *  cufflinks
 *
 *  Created by Adam Roberts on 4/5/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 *
 */

#include <assert.h>
#include <algorithm>
#include "biascorrection.h"
#include "transcripts.h"
#include "fragments.h"
#include "frequencymatrix.h"

using namespace std;

int NUM_NUCS = 4;

SeqWeightTable::SeqWeightTable(size_t window_size, double alpha)
:_observed(window_size, NUM_NUCS, alpha),
 _expected(1, NUM_NUCS, 0) 
{}

inline size_t SeqWeightTable::ctoi(char c) const
{
    switch(c)
    {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return 4;
    }
}

void SeqWeightTable::increment_expected(char c)
{
    boost::mutex::scoped_lock lock(_lock);
    size_t index = ctoi(c);
    if (index != 4)
    {
        _expected.increment(index, 1.0);
    }
}

void SeqWeightTable::increment_observed(string& seq, double normalized_mass)
{
    boost::mutex::scoped_lock lock(_lock);
    for (size_t i = 0; i < seq.length(); ++i)
    {
        size_t index = ctoi(seq[i]);
        if (index != 4)
            _observed.increment(index, normalized_mass);
    }
}

const size_t SURROUND = 10;
const size_t CENTER = 11;
const size_t WINDOW = 21;
string PADDING = "NNNNNNNNNN"; 

double SeqWeightTable::get_weight(const string& seq, size_t i) const
{
    boost::mutex::scoped_lock lock(_lock);
    double weight = 1.0;
    for (size_t j = max((size_t)0, CENTER-1-i); j < min(WINDOW, CENTER-1+seq.length()-i); ++j)
    {
        size_t index = ctoi(seq[i+j-CENTER+1]);
        if (index != 4)
            weight *= (_observed(j, index) / _expected(index));
    }
    return weight;
}

BiasBoss::BiasBoss(double alpha)
: _5_seq_bias(WINDOW, alpha),
  _3_seq_bias(WINDOW, alpha)
{}

void BiasBoss::update_expectations(const Transcript& trans)
{
    const string& seq = trans.seq();
    for (size_t i = 0; i < seq.length(); ++i)
    {
        _5_seq_bias.increment_expected(seq[i]);
        _3_seq_bias.increment_expected(seq[i]);
    }
}

void BiasBoss::update_observed(const FragMap& frag, const Transcript& trans, double normalized_mass)
{
    assert (frag.length() > WINDOW);
    
    string seq_5;
    int left_window = frag.left - (CENTER-1);
    if (left_window < 0)
    {
        seq_5 = PADDING.substr(0, -left_window);
        seq_5 += trans.seq().substr(0, WINDOW + left_window);
    }
    else
    {
        seq_5 = trans.seq().substr(left_window, WINDOW); 
    }
    _5_seq_bias.increment_observed(seq_5, normalized_mass);
    
    left_window = frag.right - CENTER;
    string seq_3 = trans.seq().substr(left_window, WINDOW);
    int overhang =left_window + WINDOW - trans.length();
    if (overhang > 0)
    {
        seq_3 += PADDING.substr(0, overhang);
    }
    _3_seq_bias.increment_observed(seq_3, normalized_mass);
}

void BiasBoss::get_transcript_bias(std::vector<double>& start_bias, std::vector<double>& end_bias, const Transcript& trans) const
{
    size_t n;
    for (size_t i = 0; i < trans.length(); ++i)
    {
        start_bias[i] = _5_seq_bias.get_weight(trans.seq(), i);
        end_bias[i] = _3_seq_bias.get_weight(trans.seq(), i);
    }
}


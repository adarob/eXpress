/*
 *  biascorrection.h
 *  express
 *
 *  Created by Adam Roberts on 4/5/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 *
 */

#include <cassert>
#include <algorithm>
#include <boost/assign.hpp>
#include "main.h"
#include "biascorrection.h"
#include "transcripts.h"
#include "fragments.h"
#include "frequencymatrix.h"
#include "sequence.h"

using namespace std;

const vector<size_t> LEN_BINS = boost::assign::list_of(791)(1265)(1707)(2433)(999999999); 
const vector<double> POS_BINS = boost::assign::list_of(0.1)(0.2)(0.3)(0.4)(0.5)(0.6)(0.7)(0.8)(0.9)(1.0);

const int SURROUND = 10;
const int CENTER = 11;
const int WINDOW = 21;
const size_t FG_ORDER = 3;
const size_t BG_ORDER = 5;

SeqWeightTable::SeqWeightTable(size_t window_size, double alpha)
:_observed(FG_ORDER, window_size, window_size, alpha),
 _expected(BG_ORDER, window_size, 1, 0, false) 
{}

void SeqWeightTable::increment_expected(const Sequence& seq)
{
    WriteLock lock(_lock);
    _expected.fast_learn(seq, 1);
}

void SeqWeightTable::normalize_expected()
{
    _expected.set_logged(true);
}

void SeqWeightTable::increment_observed(const Sequence& seq, size_t i, double normalized_mass)
{
    WriteLock lock(_lock);
    int left = (int)i - SURROUND;
    
    _observed.update(seq, left, normalized_mass);
}

double SeqWeightTable::get_weight(const Sequence& seq, size_t i) const
{
    ReadLock lock(_lock);
    int left = (int)i - SURROUND;

    return _observed.seq_prob(seq, left) - _expected.seq_prob(seq, left);
}


void SeqWeightTable::append_output(ofstream& outfile) const
{
    char buff[200];
    string header = "";
    for(int i = 0; i < WINDOW; i++)
    {
        sprintf(buff, "\t%d", i-CENTER+1);
        header += buff;
    }
    header += '\n';
    
    outfile << "\tObserved Marginal Distribution\n" << header;
    
    ReadLock lock(_lock);
    for(size_t j = 0; j < NUM_NUCS; j++)
    {
        outfile << NUCS[j] << ":\t";
        for(int i = 0; i < WINDOW; i++)
        {
            outfile << scientific << sexp(_observed.marginal_prob(i,j)) << "\t";
        }
        outfile<<endl;
    }
    
    outfile << "\tObserved Transitions\nPosition\t" << header;
    
    for(size_t j = 0; j < pow((double)NUM_NUCS, (double)FG_ORDER+1); j++)
    {
        string trans = "->" + NUCS[j & 4];
        for(size_t k = 0; k < FG_ORDER-1; ++k)
        {
            j = j >> 2;
            trans = NUCS[j & 4] + trans;
        }
        outfile << trans << '\t';
    }
    outfile << endl;
    
    for (size_t i = 0; i < WINDOW; i++)
    {
        outfile << (int)i-SURROUND << ":\t";
        for(size_t j = 0; j < pow((double)NUM_NUCS, (double)FG_ORDER+1); j++)
        {
            outfile << scientific << sexp(_observed.transition_prob(i,j>>2,j&4)) << "\t";
        }
        outfile<<endl;
    }
    
    outfile << "\tExpected Transitions\nPosition\t" << header;
    
    for(size_t j = 0; j < FG_ORDER*NUM_NUCS; j++)
    {
        string trans = "->" + NUCS[j & 4];
        for(size_t k = 0; k < FG_ORDER-1; ++k)
        {
            j = j >> 2;
            trans = NUCS[j & 4] + trans;
        }
        outfile << trans << '\t';
    }
    outfile << endl;
    
    for (size_t i = 0; i < WINDOW; i++)
    {
        outfile << (int)i-SURROUND << ":\t";
        for(size_t j = 0; j < FG_ORDER*NUM_NUCS; j++)
        {
            outfile << scientific << sexp(_expected.transition_prob(i,j>>2,j&4)) << "\t";
        }
        outfile<<endl;
    }
}

PosWeightTable::PosWeightTable(const vector<size_t>& len_bins, const vector<double>& pos_bins, double alpha)
:_observed(FrequencyMatrix(len_bins.size(), pos_bins.size(),alpha)),
 _expected(FrequencyMatrix(len_bins.size(), pos_bins.size(),0,false)),
 _len_bins(len_bins),
 _pos_bins(pos_bins)
{}

void PosWeightTable::increment_expected(size_t len, double pos)
{
    size_t l = upper_bound(len_bins().begin(),len_bins().end(), len) - len_bins().begin();
    size_t p = upper_bound(pos_bins().begin(),pos_bins().end(), pos) - pos_bins().begin();
    increment_expected(l,p);
}

void PosWeightTable::increment_expected(size_t l, size_t p)
{
    WriteLock lock(_lock);
    _expected.increment(l, p, 1);
}

void PosWeightTable::normalize_expected()
{
    _expected.set_logged(true);
}

void PosWeightTable::increment_observed(size_t len, double pos, double normalized_mass)
{
    size_t l = upper_bound(len_bins().begin(),len_bins().end(), len) - len_bins().begin();
    size_t p = upper_bound(pos_bins().begin(),pos_bins().end(), pos) - pos_bins().begin();
    increment_observed(l,p, normalized_mass);
} 

void PosWeightTable::increment_observed(size_t l, size_t p, double normalized_mass)
{
    WriteLock lock(_lock);
    _observed.increment(l,p, normalized_mass);
} 

double PosWeightTable::get_weight(size_t len, double pos) const
{
    size_t l = upper_bound(len_bins().begin(),len_bins().end(), len) - len_bins().begin();
    size_t p = upper_bound(pos_bins().begin(),pos_bins().end(), pos) - pos_bins().begin();
    ReadLock lock(_lock);
    return _observed(l,p)-_expected(l,p);
}

double PosWeightTable::get_weight(size_t l, size_t p) const
{
    ReadLock lock(_lock);
    return _observed(l,p)-_expected(l,p);
}

void PosWeightTable::append_output(ofstream& outfile) const
{
    char buff[200];
    string header = "";

    
    sprintf(buff, "\t%0.2f-%0.2f", 0.0, pos_bins()[0]);
    header += buff;
    for(size_t p = 1; p < pos_bins().size(); p++)
    {
        sprintf(buff, "\t%0.2f-%0.2f", pos_bins()[p-1], pos_bins()[p]);
        header += buff;
    }
    header += '\n';
        
    outfile << "\tObserved Position Distribution\n" << header;

    ReadLock lock(_lock);
    sprintf(buff, "%d-" SIZE_T_FMT ":\t", 0, len_bins()[0]);
    for(size_t l = 0; l < len_bins().size(); l++)
    {
        if(l)
            sprintf(buff, "" SIZE_T_FMT "-" SIZE_T_FMT ":\t", len_bins()[l-1]+1, len_bins()[l]);
        outfile << buff;
        for(size_t p = 0; p < pos_bins().size(); p++)
        {
            outfile << scientific << sexp(_observed(l,p)) << "\t";
        }
        outfile<<endl;
    }
    
    outfile << "\tBias Weights\n" << header;
    
    sprintf(buff, "%d-" SIZE_T_FMT ":\t", 0, len_bins()[0]);
    for(size_t l = 0; l < len_bins().size(); l++)
    {
        if(l)
            sprintf(buff, "" SIZE_T_FMT "-" SIZE_T_FMT ":\t", len_bins()[l-1]+1, len_bins()[l]);
        outfile << buff;
        for(size_t p = 0; p < pos_bins().size(); p++)
        {
            outfile << scientific << sexp(_observed(l,p)-_expected(l,p)) << "\t";
        }
        outfile<<endl;
    }

}

BiasBoss::BiasBoss(double alpha)
: _5_seq_bias(WINDOW, alpha),
  _3_seq_bias(WINDOW, alpha),
  _5_pos_bias(LEN_BINS, POS_BINS, alpha),
  _3_pos_bias(LEN_BINS, POS_BINS, alpha)
{}

void BiasBoss::update_expectations(const Transcript& trans)
{
    
    size_t l = upper_bound(_5_pos_bias.len_bins().begin(),_5_pos_bias.len_bins().end(), trans.length()) - _5_pos_bias.len_bins().begin();
    size_t p = 0;
    double next_bin_start = trans.length() * _5_pos_bias.pos_bins()[p];
    const Sequence& seq_fwd = trans.seq(0);
    const Sequence& seq_rev = trans.seq(1);
    for (size_t i = 0; i < trans.length(); ++i)
    {
        if (i >= next_bin_start)
        {
            next_bin_start = trans.length() * _5_pos_bias.pos_bins()[++p];
        }
        _5_pos_bias.increment_expected(l,p);
        _3_pos_bias.increment_expected(l,p);
        _5_seq_bias.increment_expected(seq_fwd);
        _3_seq_bias.increment_expected(seq_rev);
    }
}

void BiasBoss::normalize_expectations()
{
    _5_seq_bias.normalize_expected();
    _3_seq_bias.normalize_expected();
    _5_pos_bias.normalize_expected();
    _3_pos_bias.normalize_expected();
}

void BiasBoss::update_observed(const FragHit& hit, double normalized_mass)
{
    assert (hit.pair_status() != PAIRED || hit.length() > WINDOW);
    
    const Sequence& t_seq_fwd = hit.mapped_trans->seq(0);
    const Sequence& t_seq_rev = hit.mapped_trans->seq(1);

    if (hit.pair_status() != RIGHT_ONLY)
    {
        _5_seq_bias.increment_observed(t_seq_fwd, hit.left, normalized_mass);
        _5_pos_bias.increment_observed(t_seq_fwd.length(), (double)hit.left/t_seq_fwd.length(), normalized_mass);
    }
    
    if (hit.pair_status() != LEFT_ONLY)
    {
        _3_seq_bias.increment_observed(t_seq_rev, t_seq_rev.length()-hit.right, normalized_mass);
        _3_pos_bias.increment_observed(t_seq_rev.length(), (double)(hit.right-1)/t_seq_rev.length(), normalized_mass);
    }
}

double BiasBoss::get_transcript_bias(std::vector<float>& start_bias, std::vector<float>& end_bias, const Transcript& trans) const
{
    double tot_start = HUGE_VAL;
    double tot_end = HUGE_VAL;
    
    size_t l = upper_bound(_5_pos_bias.len_bins().begin(),_5_pos_bias.len_bins().end(), trans.length()) - _5_pos_bias.len_bins().begin();
    
    const Sequence& t_seq_fwd = trans.seq(0);
    const Sequence& t_seq_rev = trans.seq(1);
    
    size_t p = 0;
    double next_bin_start = trans.length() * _5_pos_bias.pos_bins()[p];
    double curr_5_pos_bias = _5_pos_bias.get_weight(l,p);
    double curr_3_pos_bias = _3_pos_bias.get_weight(l,p);
    for (size_t i = 0; i < trans.length(); ++i)
    {
        if (i >= next_bin_start)
        {
            next_bin_start = trans.length() * _5_pos_bias.pos_bins()[++p];
            curr_5_pos_bias = _5_pos_bias.get_weight(l,p);
            curr_3_pos_bias = _3_pos_bias.get_weight(l,p);
        }
        start_bias[i] = _5_seq_bias.get_weight(t_seq_fwd, i) + curr_5_pos_bias;
        end_bias[i] = _3_seq_bias.get_weight(t_seq_rev, i) + curr_3_pos_bias;
        tot_start = log_sum(tot_start, start_bias[i]);
        tot_end = log_sum(tot_start, end_bias[i]);
    }
    
    double avg_bias = (tot_start + tot_end) - (2*log((double)trans.length()));
    return avg_bias;
}

string BiasBoss::to_string() const
{  
    return "";
}

void BiasBoss::append_output(ofstream& outfile) const
{
    outfile<<">5' Sequence-Specific Bias\n";
    _5_seq_bias.append_output(outfile);
    outfile<<">3' Sequence-Specific Bias\n";
    _3_seq_bias.append_output(outfile);
    outfile<<">5' Fractional Position Bias\n";
    _5_pos_bias.append_output(outfile);
    outfile<<">3' Fractional Position Bias\n";
    _3_pos_bias.append_output(outfile);
}


/**
 *  biascorrection.h
 *  express
 *
 *  Created by Adam Roberts on 4/5/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 *
 */

#include <algorithm>
#include <boost/assign.hpp>
#include <cassert>

#include "biascorrection.h"

#include "fragments.h"
#include "frequencymatrix.h"
#include "main.h"
#include "sequence.h"
#include "targets.h"

using namespace std;

// size of bias window
const int WINDOW = 21;
// center position in window
const int CENTER = 11;
// number of bases on either side of center
const int SURROUND = 10;
// order of markov model for bias foreground
const size_t FG_ORDER = 3;
// order of markov model for bias background
const size_t BG_ORDER = 3;

SeqWeightTable::SeqWeightTable(size_t window_size, double alpha)
    : _observed(FG_ORDER, window_size, window_size, alpha),
      _expected(BG_ORDER, window_size, BG_ORDER+1, 0) {
}

void SeqWeightTable::copy_observed(const SeqWeightTable& other) {
  _observed = other._observed;
}

void SeqWeightTable::copy_expected(const SeqWeightTable& other) {
  _expected = other._expected;
}

void SeqWeightTable::increment_expected(const Sequence& seq, double mass,
                                        const vector<double>& fl_cdf) {
  _expected.fast_learn(seq, mass, fl_cdf);
}

void SeqWeightTable::normalize_expected() {
  _expected.calc_marginals();
}

void SeqWeightTable::increment_observed(const Sequence& seq, size_t i,
                                        double mass) {
  int left = (int)i - SURROUND;
  _observed.update(seq, left, mass);
}

double SeqWeightTable::get_weight(const Sequence& seq, size_t i) const {
  int left = (int)i - SURROUND;
  return _observed.seq_prob(seq, left) - _expected.seq_prob(seq, left);
}


void SeqWeightTable::append_output(ofstream& outfile) const {
  char buff[200];
  string header = "";
  for (int i = 0; i < WINDOW; i++) {
    sprintf(buff, "\t%d", i-CENTER+1);
    header += buff;
  }
  header += '\n';
    
  outfile << "\tObserved Marginal Distribution\n" << header;
    
  for (size_t j = 0; j < NUM_NUCS; j++) {
    outfile << NUCS[j] << ":\t";
    for (int i = 0; i < WINDOW; i++) {
      outfile << scientific << sexp(_observed.marginal_prob(i,j)) << "\t";
    }
    outfile<<endl;
  }
    
  outfile << "\tObserved Conditional Probabilities\nPosition\t";
    
  for (size_t j = 0; j < pow((double)NUM_NUCS, (double)FG_ORDER+1); j++) {
    string s = "->";
    s += NUCS[j & 3];
    size_t cond = j >> 2;
    for (size_t k = 0; k < FG_ORDER; ++k) {
      s = NUCS[cond & 3] + s;
      cond = cond >> 2;
    }
    outfile << s << '\t';
  }
  outfile << endl;
    
  for (int i = 0; i < WINDOW; i++) {
    outfile << i-SURROUND << ":\t";
    for (size_t j = 0; j < pow((double)NUM_NUCS, (double)FG_ORDER+1); j++) {
      outfile << scientific << sexp(_observed.transition_prob(i,j>>2,j&3))
              << "\t";
    }
    outfile<<endl;
  }
    
  outfile << "\tBackground Conditional Probabilities\nPosition\t";
    
  for (size_t j = 0; j < pow((double)NUM_NUCS, (double)BG_ORDER+1); j++) {
    string s = "->";
    s += NUCS[j & 3];
    size_t cond = j >> 2;
    for (size_t k = 0; k < BG_ORDER; ++k) {
      s = NUCS[cond & 3] + s;
      cond = cond >> 2;
    }
    outfile << s << '\t';
  }
  outfile << endl;
  
  for (size_t i = 0; i < BG_ORDER+1; i++) {
    outfile << i << ":\t";
    for (size_t j = 0; j <  pow((double)NUM_NUCS, (double)BG_ORDER+1); j++) {
      outfile << scientific << sexp(_expected.transition_prob(i,j>>2,j&3))
              << "\t";
    }
    outfile<<endl;
  }
}

BiasBoss::BiasBoss(double alpha)
    : _5_seq_bias(WINDOW, alpha),
      _3_seq_bias(WINDOW, alpha){
}

void BiasBoss::copy_observations(const BiasBoss& other) {
  _5_seq_bias.copy_observed(other._5_seq_bias);
  _3_seq_bias.copy_observed(other._3_seq_bias);
}

void BiasBoss::copy_expectations(const BiasBoss& other) {
  _5_seq_bias.copy_expected(other._5_seq_bias);
  _3_seq_bias.copy_expected(other._3_seq_bias);
}

void BiasBoss::update_expectations(const Target& targ, double mass,
                                   const vector<double>& fl_cdf) {
  if (mass == LOG_0) {
    return;
  }
  
  const Sequence& seq_fwd = targ.seq(0);
  const Sequence& seq_rev = targ.seq(1);

  _5_seq_bias.increment_expected(seq_fwd, mass, fl_cdf);
  _3_seq_bias.increment_expected(seq_rev, mass, fl_cdf);
}

void BiasBoss::normalize_expectations() {
  _5_seq_bias.normalize_expected();
  _3_seq_bias.normalize_expected();
}

void BiasBoss::update_observed(const FragHit& hit, double normalized_mass)
{
  assert (hit.pair_status() != PAIRED || hit.length() > WINDOW);
    
  const Sequence& t_seq_fwd = hit.targ->seq(0);
  const Sequence& t_seq_rev = hit.targ->seq(1);

  if (hit.pair_status() != RIGHT_ONLY) {
    _5_seq_bias.increment_observed(t_seq_fwd, hit.left, normalized_mass);
  }
  if (hit.pair_status() != LEFT_ONLY) {
    _3_seq_bias.increment_observed(t_seq_rev, t_seq_rev.length()-hit.right,
                                   normalized_mass);
  }
}

double BiasBoss::get_target_bias(std::vector<float>& start_bias,
                                 std::vector<float>& end_bias,
                                 const Target& targ) const {
  double tot_start = LOG_0;
  double tot_end = LOG_0;
    
  const Sequence& t_seq_fwd = targ.seq(0);
  const Sequence& t_seq_rev = targ.seq(1);

  for (size_t i = 0; i < targ.length(); ++i) {
    start_bias[i] = _5_seq_bias.get_weight(t_seq_fwd, i);
    end_bias[targ.length()-i-1] = _3_seq_bias.get_weight(t_seq_rev, i);
    tot_start = log_add(tot_start, start_bias[i]);
    tot_end = log_add(tot_start, end_bias[i]);
  }
    
  double avg_bias = (tot_start + tot_end) - (2*log((double)targ.length()));
  assert(!isnan(avg_bias));
  return avg_bias;
}

void BiasBoss::append_output(ofstream& outfile) const {
  outfile << ">5' Sequence-Specific Bias\n";
  _5_seq_bias.append_output(outfile);
  outfile << ">3' Sequence-Specific Bias\n";
  _3_seq_bias.append_output(outfile);
}

//
//  sequence.cpp
//  express
//
//  Created by Adam Roberts on 1/24/12.
//  Copyright 2012 Adam Roberts. All rights reserved.
//

#include "sequence.h"
#include <cassert>
#include <boost/math/distributions/normal.hpp>

using namespace std;
using namespace boost::math;

string Sequence::serialize() {
  vector<char> seq;
  for (size_t i = 0; i < length(); i++) {
    if (i/4 == seq.size()) {
      seq.push_back(0);
    }
    seq.back() += operator[](i) << (2 * (i % 4));
  }
  
  return string(seq.begin(), seq.end());
}

SequenceFwd::SequenceFwd():  _ref_seq(NULL), _prob(0), _len(0) {}

SequenceFwd::SequenceFwd(const std::string& seq, bool rev, bool prob)
    : _prob(prob), _len(seq.length()) {
  if (prob) {
    _est_seq = FrequencyMatrix<float>(seq.length(), NUM_NUCS, 0.001);
    _obs_seq = FrequencyMatrix<float>(seq.length(), NUM_NUCS, LOG_0);
    _exp_seq = FrequencyMatrix<float>(seq.length(), NUM_NUCS, LOG_0);
  }
  set(seq, rev);
}

SequenceFwd::SequenceFwd(const SequenceFwd& other)
    : _obs_seq(other._obs_seq), _exp_seq(other._exp_seq),
      _prob(other._prob), _len(other.length()) {
  if (other._ref_seq) {
    char* ref_seq = new char[_len];
    std::copy(other._ref_seq.get(), other._ref_seq.get() + _len, ref_seq);
    _ref_seq.reset(ref_seq);
  }
}

SequenceFwd& SequenceFwd::operator=(const SequenceFwd& other) {
  if (other._ref_seq) {
    _len = other.length();
    char* ref_seq = new char[_len];
    std::copy(other._ref_seq.get(), other._ref_seq.get() + _len, ref_seq);
    _ref_seq.reset(ref_seq);
    _obs_seq = other._obs_seq;
    _exp_seq = other._exp_seq;
    _prob = other._prob;
  }
  return *this;
}

void SequenceFwd::set(const std::string& seq, bool rev) {
  char* ref_seq = new char[seq.length()];
  for (size_t i = 0; i < seq.length(); i++) {
    ref_seq[i] = (rev) ? complement(ctoi(seq[seq.length()-1-i])) : ctoi(seq[i]);
    if (_prob) {
      _est_seq.increment(i, ref_seq[i], log((float)2));
    }
  }
  _ref_seq.reset(ref_seq);
  _len = seq.length();
}

size_t SequenceFwd::operator[](const size_t index) const {
  assert(index < _len);
  if (_prob) {
    return _est_seq.argmax(index);
  }
  return _ref_seq[index];
}

size_t SequenceFwd::get_ref(const size_t index) const {
  assert(index < _len);
  //assert(_ref_seq[index] == operator[](index));
  return _ref_seq[index];
}

float SequenceFwd::get_prob(const size_t index, const size_t nuc) const {
  assert(_prob);
  return _est_seq(index, nuc);
}

float SequenceFwd::get_obs(const size_t index, const size_t nuc) const {
  assert(index < _len);
  return _obs_seq(index,nuc, false);
}

float SequenceFwd::get_exp(const size_t index, const size_t nuc) const {
  assert(index < _len);
  return _exp_seq(index,nuc, false);
}

void SequenceFwd::update_est(const size_t index, const size_t nuc, float mass) {
  assert(_prob);
  _est_seq.increment(index, nuc, mass);
}

void SequenceFwd::update_obs(const size_t index, const size_t nuc, float mass) {
  assert(_prob);
  _obs_seq.increment(index, nuc, mass);
}

void SequenceFwd::update_exp(const size_t index, const size_t nuc, float mass) {
  assert(_prob);
  _exp_seq.increment(index, nuc, mass);
}

void SequenceFwd::calc_p_vals(vector<double>& p_vals) const {
  p_vals = vector<double>(_len, 1.0);
  for (size_t i = 0; i < _len; ++i) {
    double N = sexp(_obs_seq.sum(i));
    if (N==0) {
      continue;
    }
    size_t ref_nuc = get_ref(i);
    double max_obs = 0;
    for (size_t nuc = 0; nuc < NUM_NUCS; ++nuc) {
      if (nuc == ref_nuc) {
        continue;
      }

      double obs_n = sexp(_obs_seq(i,nuc,false));
      max_obs = max(max_obs,obs_n);
    }

    double p_val = 0;

    for (size_t nuc = 0; nuc < NUM_NUCS; ++nuc) {
      if (nuc == ref_nuc) {
        continue;
      }

      double exp_p = sexp(_exp_seq(i, nuc));
      normal norm(N*exp_p, sqrt(N*exp_p*(1-exp_p)));
      p_val += log(cdf(norm, max_obs));
    }
    p_vals[i] -= sexp(p_val);
  }
}

void SequenceRev::calc_p_vals(vector<double>& p_vals) const
{
  vector<double> temp;
  _seq->calc_p_vals(temp);
  p_vals = vector<double>(length());
  for(size_t i = 0; i < length(); i++) {
    p_vals[i] = temp[length()-i-1];
  }
}

/*
void SequenceFwd::calc_p_vals(vector<double>& p_vals) const
{
    p_vals = vector<double>(_len, 1);
    for (size_t i = 0; i < _len; ++i)
    {
        double N = sexp(_obs_seq.sum(i));

        if (N==0)
            continue;

        size_t ref_nuc = get_ref(i);
        double p_val = LOG_0;

        vector<double> cdfs(4);
        for (size_t nuc = 0; nuc < NUM_NUCS; ++nuc)
        {
            if (nuc == ref_nuc)
                continue;

            double p = sexp(_exp_seq(i,nuc));
            double obs_n = sexp(_obs_seq(i,nuc,false));
            normal norm (N*p, sqrt(N*p*(1-p)));
            cdfs[nuc] = cdf(norm, obs_n);
        }

        for (size_t nuc1 = 0; nuc1 < NUM_NUCS; ++nuc1)
        {
            if (nuc1 == ref_nuc)
                continue;

            double term = log(1-cdfs[nuc1]);

            for(size_t nuc2 = 0; nuc2 < NUM_NUCS; ++nuc2)
            {
                if (nuc2 == nuc1 || nuc2 == ref_nuc)
                    continue;
                term += log(cdfs[nuc2]);
            }
            p_val = log_add(p_val, term);
        }

        p_vals[i] = sexp(p_val);
    }
}
*/

/*
void SequenceFwd::calc_p_vals(vector<double>& p_vals) const
{
    p_vals = vector<double>(_len, 1.0);
    boost::math::chi_squared_distribution<double> chisq(3);
    double obs_n,exp_n;
    for (size_t i = 0; i < _len; ++i)
    {
        if (_obs_seq.sum(i)==LOG_0)
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
 */

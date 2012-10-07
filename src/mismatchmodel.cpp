//
//  mismatchmodel.cpp
//  express
//
//  Created by Adam Roberts on 4/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "main.h"
#include "mismatchmodel.h"
#include "targets.h"
#include "fragments.h"
#include "sequence.h"
#include <iostream>
#include <fstream>

using namespace std;

MismatchTable::MismatchTable(double alpha)
    : _first_read_mm(MAX_READ_LEN, FrequencyMatrix<double>(16, 4, alpha)),
      _second_read_mm(MAX_READ_LEN, FrequencyMatrix<double>(16, 4, alpha)),
      _insert_params(1, MAX_READ_LEN, alpha),
      _delete_params(1, MAX_READ_LEN, alpha),
      _max_len(0),
      _active(false){
}

void MismatchTable::get_indices(const FragHit& f,
                           vector<char>& left_indices,
                           vector<char>& right_indices) const {

  const Target& targ = *f.target();
  const Sequence& t_seq_fwd = targ.seq(0);
  const Sequence& t_seq_rev = targ.seq(1);
  
  if (f.left_read()) {
    const ReadHit& read_l = *f.left_read();
    
    left_indices = vector<char>(read_l.seq.length(), 255);
    
    size_t i = 0;  // read index
    size_t j = read_l.left;  // genomic index
    
    vector<Indel>::const_iterator ins = read_l.inserts.begin();
    vector<Indel>::const_iterator del = read_l.deletes.begin();
    
    while (i < read_l.seq.length()) {
      if (del != read_l.deletes.end() && del->pos == i) {
        j += del->len;
        del++;
      } else if (ins != read_l.inserts.end() && ins->pos == i) {
        i += ins->len;
        ins++;
      } else {
        size_t cur = read_l.seq[i];
        size_t prev = (i) ? (read_l.seq[i-1] << 2) : 0;
        size_t ref = t_seq_fwd[j];
        size_t index = prev + ref;
        left_indices[i] = (index << 2) + cur;
        
        i++;
        j++;
      }
    }
  }
  
  if (f.right_read()) {
    const ReadHit& read_r = *f.right_read();
    
    right_indices = vector<char>(read_r.seq.length(), 255);

    size_t r_len = read_r.seq.length();
    size_t i = 0;
    size_t j = targ.length() - read_r.right;
        
    vector<Indel>::const_iterator ins = read_r.inserts.end()-1;
    vector<Indel>::const_iterator del = read_r.deletes.end()-1;
    
    while (i < r_len) {
      if (del != read_r.deletes.begin() - 1 && del->pos == r_len - i) {
        j += del->len;
        del--;
      } else if (ins != read_r.inserts.begin() - 1
                 && ins->pos + ins->len == r_len - i) {
        i += ins->len;
        ins--;
      } else {
        size_t cur = read_r.seq[i];
        size_t prev = (i) ? (read_r.seq[i-1] << 2) : 0;
        size_t ref = t_seq_rev[j];
        size_t index = prev + ref;;
        right_indices[i] = (index << 2) + cur;
        
        i++;
        j++;
      }
    }
  }
}

double MismatchTable::log_likelihood(const FragHit& f) const {
  if (!_active) {
    return 0;
  }
  
  const Target& targ = *f.target();
  const Sequence& t_seq_fwd = targ.seq(0);
  const Sequence& t_seq_rev = targ.seq(1);

  double ll = 0;

  if (f.left_read()) {
    const ReadHit& read_l = *f.left_read();
    const vector<FrequencyMatrix<double> >& left_mm = (read_l.first) ?
                                                      _first_read_mm :
                                                      _second_read_mm;
    size_t i = 0;  // read index
    size_t j = read_l.left;  // genomic index

    size_t del_len = 0;
    vector<Indel>::const_iterator ins = read_l.inserts.begin();
    vector<Indel>::const_iterator del = read_l.deletes.begin();
    
    while (i < read_l.seq.length()) {
      if (del != read_l.deletes.end() && del->pos == i) {
        ll += _delete_params(del->len);
        j += del->len;
        del_len = del->len;
        del++;
      } else if (ins != read_l.inserts.end() && ins->pos == i) {
        ll += _insert_params(ins->len);
        i += ins->len;
        if (ins->len > del_len) {
          ll += (ins->len-del_len)*_delete_params(0);  //FIXME: Assymetric
        }
        del_len = 0;
        ins++;
      } else {
        if (del_len > 0) {
          ll += _delete_params(0);
        }
        del_len = 0;
        ll += _insert_params(0);

        size_t cur = read_l.seq[i];
        size_t prev = (i) ? (read_l.seq[i-1] << 2) : 0;

        if (t_seq_fwd.prob()) {
          double trans_prob = LOG_0;
          for (size_t nuc = 0; nuc < NUM_NUCS; nuc++) {
            size_t index = ((prev + nuc) << 2) + cur;
            
            trans_prob = log_add(trans_prob, t_seq_fwd.get_prob(j, nuc) +
                                             left_mm[i](index));
          }
          ll += trans_prob;
        } else {
          size_t ref = t_seq_fwd[j];
          size_t index = prev + ref;
          ll += left_mm[i](index, cur);
        }
        i++;
        j++;
      }
    }
  }
  
  if (f.right_read()) {
    const ReadHit& read_r = *f.right_read();
    
    const vector<FrequencyMatrix<double> >& right_mm = (read_r.first) ?
                                                        _first_read_mm :
                                                        _second_read_mm;
    
    size_t r_len = read_r.seq.length();
    size_t i = 0;
    size_t j = targ.length() - read_r.right;

    size_t del_len = 0;

    vector<Indel>::const_iterator ins = read_r.inserts.end()-1;
    vector<Indel>::const_iterator del = read_r.deletes.end()-1;

    while (i < r_len) {
      if (del != read_r.deletes.begin() - 1 && del->pos == r_len - i) {
        ll += _delete_params(del->len);
        j += del->len;
        del_len = del->len;
        del--;
      } else if (ins != read_r.inserts.begin() - 1
                 && ins->pos + ins->len == r_len - i) {
        ll += _insert_params(ins->len);
        if (ins->len > del_len) {
          ll += (ins->len-del_len) * _delete_params(0);
        }
        del_len = 0;

        i += ins->len;
        ins--;
      } else {
        if (del_len > 0) {
          ll += _delete_params(0);
        }
        del_len = 0;
        ll += _insert_params(0);

        size_t cur = read_r.seq[i];
        size_t prev = (i) ? (read_r.seq[i-1] << 2) : 0;

        if (t_seq_rev.prob()) {
          double trans_prob = LOG_0;
          for (size_t nuc = 0; nuc < NUM_NUCS; nuc++) {
            size_t index = ((prev + nuc) << 2) + cur;
            trans_prob = log_add(trans_prob, t_seq_rev.get_prob(j, nuc) +
                                             right_mm[i](index));
          }
          ll += trans_prob;
        } else {
          size_t ref = t_seq_rev[j];
          size_t index = prev + ref;;
          ll += right_mm[i](index, cur);
        }
        i++;
        j++;
      }
    }
  }
  
  assert(!(isnan(ll)||isinf(ll)));
  return ll;
}

void MismatchTable::update(const FragHit& f, double p, double mass) {
  Target& targ = *f.target();
  Sequence& t_seq_fwd = targ.seq(0);
  Sequence& t_seq_rev = targ.seq(1);

  if (f.left_read()) {
    const ReadHit& read_l = *f.left_read();
    vector<FrequencyMatrix<double> >& left_mm = (read_l.first) ?
                                                _first_read_mm :
                                                _second_read_mm;
    size_t i = 0;  // read index
    size_t j = read_l.left;  // genomic index

    size_t del_len = 0;
    vector<Indel>::const_iterator ins = read_l.inserts.begin();
    vector<Indel>::const_iterator del = read_l.deletes.begin();

    vector<double> joint_probs(NUM_NUCS);

    assert(targ.length() >= f.right());
    while (i < read_l.seq.length()) {
      if (del != read_l.deletes.end() && del->pos == i) {
        _delete_params.increment(del->len, mass);
        j += del->len;
        del_len = del->len;
        del++;
      } else if (ins != read_l.inserts.end() && ins->pos == i) {
        _insert_params.increment(ins->len, mass);
        i += ins->len;
        if (ins->len > del_len) {
          _delete_params.increment(0, mass);  //FIXME: Assymetric
        }
        del_len = 0;
        ins++;
      } else {
        if (del_len > 0) {
           _delete_params.increment(0, mass);
        }
        del_len = 0;
        _insert_params.increment(0, mass);

        size_t cur = read_l.seq[i];
        size_t prev = (i) ? (read_l.seq[i-1] << 2) : 0;
        // Update the seq parameters only after burn-in (active)
        if (t_seq_fwd.prob() && _active) {
          
          t_seq_fwd.update_obs(j, cur, p);

          double Z = LOG_0;

          size_t ref_index = (i) ? (t_seq_fwd.get_ref(j-1)<<2) +
                                    t_seq_fwd.get_ref(j) : t_seq_fwd.get_ref(j);
          for (size_t nuc = 0; nuc < NUM_NUCS; nuc++) {
            // Update expected
            t_seq_fwd.update_exp(j, nuc, p+left_mm[i](ref_index, nuc));

            // Update posterior
            size_t index = prev + nuc;
            joint_probs[nuc] = t_seq_fwd.get_prob(j, nuc) +
                               left_mm[i](index, cur);
            Z = log_add(Z, joint_probs[nuc]);
          }

          for (size_t nuc = 0; !left_mm[i].is_fixed() && nuc < NUM_NUCS; nuc++) {
            size_t index = prev + nuc;
            left_mm[i].increment(index, cur,
                                 mass + p + t_seq_fwd.get_prob(j, nuc));
          }

          for (size_t nuc=0; nuc < NUM_NUCS; nuc++) {
            t_seq_fwd.update_est(j, nuc, p + joint_probs[nuc] - Z);
          }
        } else {
          size_t ref = t_seq_fwd[j];
          size_t index = prev + ref;
          left_mm[i].increment(index, cur, mass + p);
        }

        i++;
        j++;
      }
    }
    _max_len = max(_max_len, read_l.seq.length());
  }
  
  if (f.right_read()) {
    const ReadHit& read_r = *f.right_read();
    vector<FrequencyMatrix<double> >& right_mm = (read_r.first) ?
                                                  _first_read_mm :
                                                  _second_read_mm;
    
    size_t r_len = read_r.seq.length();
    size_t i = 0;
    size_t j = targ.length() - read_r.right;

    size_t del_len = 0;
    vector<Indel>::const_iterator ins = read_r.inserts.end() - 1;
    vector<Indel>::const_iterator del = read_r.deletes.end() - 1;

    vector<double> joint_probs(NUM_NUCS);

    while (i < r_len) {
      if (del != read_r.deletes.begin()-1 && del->pos == r_len-i ) {
         _delete_params.increment(del->len, mass);
         j += del->len;
         del_len = del->len;
         del--;
      } else if (ins != read_r.inserts.begin() - 1 &&
                 ins->pos + ins->len == r_len-i) {
        _insert_params.increment(ins->len, mass);
        if (ins->len > del_len) {
          _delete_params.increment(0, mass);
        }
        del_len = 0;

        i += ins->len;
        ins--;
      } else {
        if (del_len > 0) {
          _delete_params.increment(0, mass);
        }
        del_len = 0;
        _insert_params.increment(0, mass);

        size_t cur = read_r.seq[i];
        size_t prev = (i) ? (read_r.seq[i-1] << 2) : 0;

        if (t_seq_rev.prob() && _active) {          
          t_seq_rev.update_obs(j, cur, p);

          double Z = LOG_0;

          size_t ref_index = (i) ? (t_seq_rev.get_ref(j-1)<<2) +
                                   t_seq_rev.get_ref(j) : t_seq_rev.get_ref(j);
          for (size_t nuc = 0; nuc < NUM_NUCS; nuc++) {
            // Update expected
            t_seq_rev.update_exp(j, nuc, p+right_mm[i](ref_index, nuc));

            // Update posterior
            size_t index = prev + nuc;
            joint_probs[nuc] = t_seq_rev.get_prob(j, nuc) +
                               right_mm[i](index, cur);
            Z = log_add(Z, joint_probs[nuc]);
          }

          for (size_t nuc = 0; !right_mm[i].is_fixed() && nuc < NUM_NUCS; nuc++) {
            size_t index = prev + nuc;
            right_mm[i].increment(index, cur, mass+p+t_seq_rev.get_prob(j, nuc));
          }

          for (size_t nuc=0; nuc < NUM_NUCS; nuc++) {
            t_seq_rev.update_est(j, nuc, p + joint_probs[nuc] - Z);
          }
        } else {
          size_t ref = t_seq_rev[j];
          size_t index = prev + ref;
          right_mm[i].increment(index, mass+p);
        }

        i++;
        j++;
      }
      _max_len = max(_max_len, read_r.seq.length());
    }
  }
}

void MismatchTable::fix() {
  for (size_t i = 0; i < MAX_READ_LEN; i++) {
    _first_read_mm[i].fix();
    _second_read_mm[i].fix();
  }
  _insert_params.fix();
  _delete_params.fix();
}

void MismatchTable::append_output(ofstream& outfile) const {
  string col_header =  "\t";
  for (size_t i = 0; i < 64; i++) {
    col_header += NUCS[i>>4];
    col_header += NUCS[i>>2 & 3];
    col_header += "->*";
    col_header += NUCS[i & 3];
    col_header += '\t';
  }
  col_header[col_header.length()-1] = '\n';

  outfile << ">First Read Mismatch\n" << col_header;
  for (size_t k = 0; k < _max_len; k++) {
    outfile << k << ":\t";
    for (size_t i = 0; i < 16; i++) {
      for (size_t j = 0; j < 4; j++) {
        if (k || i < 4) {
          outfile << scientific << sexp(_first_read_mm[k](i,j))<<"\t";
        } else {
          outfile << scientific << 0.0 << "\t";
        }
      }
    }
    outfile<<endl;
  }

  outfile<<">Second Read Mismatch\n" << col_header;
  for (size_t k = 0; k < _max_len; k++) {
    outfile << k << ":\t";
    for (size_t i = 0; i < 16; i++) {
      for (size_t j = 0; j < 4; j++) {
        if (k || i < 4) {
          outfile << scientific << sexp(_second_read_mm[k](i,j))<<"\t";
        } else {
          outfile << scientific << 0.0 << "\t";
        }
      }
    }
    outfile<<endl;
  }
}

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

double MismatchTable::log_likelihood(const FragHit& f) const {
  if (!_active) {
    return 0;
  }
    
  const Target& targ = *f.targ;
  const Sequence& t_seq_fwd = targ.seq(0);
  const Sequence& t_seq_rev = targ.seq(1);

  double ll = 0;
    
  const vector<FrequencyMatrix<double> >& left_mm = (f.left_first) ?
                                                    _first_read_mm :
                                                    _second_read_mm;
  const vector<FrequencyMatrix<double> >& right_mm = (!f.left_first) ?
                                                      _first_read_mm :
                                                      _second_read_mm;
   
  size_t cur;
  size_t ref;
  size_t prev;
  size_t index;
  
  size_t i = 0;  // read index
  size_t j = f.left;  // genomic index
    
  size_t del_len = 0;
  vector<Indel>::const_iterator ins = f.inserts_l.begin();
  vector<Indel>::const_iterator del = f.deletes_l.begin();
    
  while (i < f.seq_l.length()) {
    if (del != f.deletes_l.end() && del->pos == i) {
      ll += _delete_params(del->len);
      j += del->len;
      del_len = del->len;
      del++;
    } else if (ins != f.inserts_l.end() && ins->pos == i) {
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
            
      cur = f.seq_l[i];
      prev = (i) ? (f.seq_l[i-1] << 2) : 0;
            
      if (t_seq_fwd.prob()) {
        double trans_prob = LOG_0;
        for (size_t nuc = 0; nuc < NUM_NUCS; nuc++) {
          index = prev + nuc;
          trans_prob = log_add(trans_prob, t_seq_fwd.get_prob(j, nuc) +
                                           left_mm[i](index, cur));
        }
        ll += trans_prob;
      } else {
        ref = t_seq_fwd[j];
        index = prev + ref;
        ll += left_mm[i](index, cur);
      }
      i++;
      j++;
    }
  }
    
  size_t r_len = f.seq_r.length();
  i = 0;
  j = targ.length()-f.right;
  ins = f.inserts_r.end()-1;
  del = f.deletes_r.end()-1;
  
  while (i < r_len) {
    if (del != f.deletes_r.begin()-1 && del->pos == r_len-i) {
      ll += _delete_params(del->len);
      j += del->len;
      del_len = del->len;
      del--;
    } else if (ins != f.inserts_r.begin()-1 && ins->pos + ins->len == r_len-i) {
      ll += _insert_params(ins->len);
      if (ins->len > del_len) {
                ll += (ins->len-del_len)*_delete_params(0);
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

      cur = f.seq_r[i];
      prev = (i) ? (f.seq_r[i-1] << 2) : 0;
            
      if (t_seq_rev.prob()) {
        double trans_prob = LOG_0;
        for (size_t nuc = 0; nuc < NUM_NUCS; nuc++) {
          index = prev + nuc;
          trans_prob = log_add(trans_prob, t_seq_rev.get_prob(j, nuc) +
                                           right_mm[i](index, cur));
        }
        ll += trans_prob;
      } else {
        ref = t_seq_rev[j];
        index = prev + ref;
        ll += right_mm[i](index, cur);
      }

      i++;
      j++;
    }
  }
    
  assert(!(isnan(ll)||isinf(ll)));
  return ll;
}


void MismatchTable::update(const FragHit& f, double p, double mass) {
  Target& targ = *f.targ;
  Sequence& t_seq_fwd = targ.seq(0);
  Sequence& t_seq_rev = targ.seq(1);
    
  vector<FrequencyMatrix<double> >& left_mm = (f.left_first) ? _first_read_mm :
                                                               _second_read_mm;
  vector<FrequencyMatrix<double> >& right_mm = (!f.left_first) ?
                                                _first_read_mm :
                                                _second_read_mm;
    
  size_t cur;
  size_t prev;
  size_t index;
    
  size_t i = 0;  // read index
  size_t j = f.left;  // genomic index
    
  size_t del_len = 0;
  vector<Indel>::const_iterator ins = f.inserts_l.begin();
  vector<Indel>::const_iterator del = f.deletes_l.begin();
    
  vector<double> joint_probs(NUM_NUCS);
    
  while (i < f.seq_l.length()) {
    if (del != f.deletes_l.end() && del->pos == i) {
      _delete_params.increment(del->len, mass);
      j += del->len;
      del_len = del->len;
      del++;
    } else if (ins != f.inserts_l.end() && ins->pos == i) {
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
            
      cur = f.seq_l[i];
      prev = (i) ? (f.seq_l[i-1] << 2) : 0;
     
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
          index = prev + nuc;
          joint_probs[nuc] = t_seq_fwd.get_prob(j, nuc) +
                             left_mm[i](index, cur);
          Z = log_add(Z, joint_probs[nuc]);
        }
                
        for (size_t nuc = 0; !left_mm[i].is_fixed() && nuc < NUM_NUCS; nuc++) {
          index = prev + nuc;
          left_mm[i].increment(index, cur, mass+p+t_seq_fwd.get_prob(j, nuc));
        }
                
        for (size_t nuc=0; nuc < NUM_NUCS; nuc++) {
          t_seq_fwd.update_est(j, nuc, p + joint_probs[nuc] - Z);
        }
      } else {
        index = prev + t_seq_fwd[j];
        left_mm[i].increment(index, cur, mass+p);
      }
            
      i++;
      j++;
    }
  }

  size_t r_len = f.seq_r.length();
  i = 0;
  j = targ.length()-f.right;
  ins = f.inserts_r.end()-1;
  del = f.deletes_r.end()-1;

  while (i < r_len) {
    if (del != f.deletes_r.begin()-1 && del->pos == r_len-i ) {
            _delete_params.increment(del->len, mass);
            j += del->len;
            del_len = del->len;
            del--;
    } else if (ins != f.inserts_r.begin()-1 && ins->pos + ins->len == r_len-i) {
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
           
      cur = f.seq_r[i];
      prev = (i) ? (f.seq_r[i-1] << 2) : 0;
            
      if (t_seq_rev.prob() && _active) {
        t_seq_rev.update_obs(j, cur, p);
               
        double Z = LOG_0;
                
        size_t ref_index = (i) ? (t_seq_rev.get_ref(j-1)<<2) +
                                 t_seq_rev.get_ref(j) : t_seq_rev.get_ref(j);
        for (size_t nuc = 0; nuc < NUM_NUCS; nuc++) {
          // Update expected
          t_seq_rev.update_exp(j, nuc, p+right_mm[i](ref_index, nuc));
                   
          // Update posterior
          index = prev + nuc;
          joint_probs[nuc] = t_seq_rev.get_prob(j, nuc) +
                             right_mm[i](index, cur);
          Z = log_add(Z, joint_probs[nuc]);
        }
                
        for (size_t nuc = 0; !right_mm[i].is_fixed() && nuc < NUM_NUCS; nuc++) {
          index = prev + nuc;
          right_mm[i].increment(index, cur, mass+p+t_seq_rev.get_prob(j, nuc));
        }
                
        for (size_t nuc=0; nuc < NUM_NUCS; nuc++) {
          t_seq_rev.update_est(j, nuc, p + joint_probs[nuc] - Z);
        }
      } else {
        index = prev + t_seq_rev[j];
        right_mm[i].increment(index, cur, mass+p);
      }
            
      i++;
      j++;
    }
  }
    
  _max_len = max(_max_len, max(f.seq_l.length(), f.seq_r.length()));
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

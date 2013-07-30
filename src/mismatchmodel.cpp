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
    : _first_read_mm(max_read_len, FrequencyMatrix<double>(16, 4, alpha)),
      _second_read_mm(max_read_len, FrequencyMatrix<double>(16, 4, alpha)),
      _insert_params(1, max_indel_size + 1, 0),
      _delete_params(1, max_indel_size + 1, 0),
      _max_len(0),
      _active(false){
  // Set indel priors
  double no_indel_p = 0.99;
  double pm = no_indel_p;
  for(size_t i = 0 ; i <= max_indel_size; ++i) {
    _insert_params.increment(i, log(alpha * pm));
    _delete_params.increment(i, log(alpha * pm));
    pm *= (1 - no_indel_p);
  }
  assert(approx_eq(sexp(_insert_params.sum(0)), alpha));
  assert(approx_eq(sexp(_delete_params.sum(0)), alpha));
}

MismatchTable::MismatchTable(string param_file_name)
    : _first_read_mm(max_read_len, FrequencyMatrix<double>(16, 4, 0)),
      _second_read_mm(max_read_len, FrequencyMatrix<double>(16, 4, 0)),
      _insert_params(1, max_indel_size + 1, 0),
      _delete_params(1, max_indel_size + 1, 0),
      _max_len(0),
      _active(true){
  ifstream infile (param_file_name.c_str());
  const size_t BUFF_SIZE = 99999;
  char line_buff[BUFF_SIZE];
  if (!infile.is_open()) {
    logger.severe("Unable to open parameter file '%s'.",
                  param_file_name.c_str());
  }

  infile.getline (line_buff, BUFF_SIZE, '\n');
  size_t pos = 0;
  
  while (infile.good()) {
    infile.getline (line_buff, BUFF_SIZE, '\n');
    if (!strcmp(line_buff, ">First Read Mismatch")) {
      break;
    }
  }
  
  infile.getline (line_buff, BUFF_SIZE, '\n');
  while (infile.good()) {
    infile.getline (line_buff, BUFF_SIZE, '\n');
    if (!strcmp(line_buff, ">Second Read Mismatch")) {
      break;
    }
    
    if (pos >= max_read_len) {
      pos++;
      continue;
    }
    
    _first_read_mm[pos] = FrequencyMatrix<double>(16, 4, 0);
    char *p = strtok(line_buff, "\t");
    for (size_t i = 0; i < 16; ++i) {
      for(size_t j = 0; j < 4; ++j) {
        p = strtok(NULL, "\t");
        _first_read_mm[pos].increment(i, j, log(strtod(p,NULL)));
      }
    }
    pos++;
  }
  _max_len = min(pos, max_read_len);
  if (pos >= max_read_len) {
    logger.warn("First read error distribution of %d bases in '%s' truncated "
                "after %d bases.",
                pos-1, param_file_name.c_str(), max_read_len);
  }
  
  pos = 0;
  infile.getline (line_buff, BUFF_SIZE, '\n');
  
  while ( infile.good() ) {
    infile.getline (line_buff, BUFF_SIZE, '\n');
    
    if (!strncmp(line_buff, ">Insertion Length", 17)) {
      break;
    }
    
    if (pos >= max_read_len) {
      pos++;
      continue;
    }
    
    _second_read_mm[pos] = FrequencyMatrix<double>(16, 4, 0);
    char *p = strtok(line_buff, "\t");
    for (size_t i = 0; i < 16; ++i) {
      for(size_t j = 0; j < 4; ++j) {
        p = strtok(NULL, "\t");
        _second_read_mm[pos].increment(i, j, log(strtod(p,NULL)));
      }
    }
    pos++;
  }
  _max_len = max(_max_len, min(pos, max_read_len));
  if (pos >= max_read_len) {
    logger.warn("Second read error distribution of %d bases in '%s' truncated "
                "after %d bases.",
                pos-1, param_file_name.c_str(), max_read_len);
  }
  
  infile.getline (line_buff, BUFF_SIZE, '\n');
  char *p = strtok(line_buff, "\t");
  size_t k = 0;
  do {
    if (k > max_indel_size) {
      logger.warn("Paramater file '%s' insertion distribution is being "
                  "truncated at max indel length of %d.",
                  param_file_name.c_str(), max_indel_size);
      break;
    }
    _insert_params.increment(k, log(strtod(p,NULL)));
    p = strtok(NULL, "\t");
    k++;

  } while (p);

  infile.getline (line_buff, BUFF_SIZE, '\n');
  infile.getline (line_buff, BUFF_SIZE, '\n');
  p = strtok(line_buff, "\t");
  k = 0;
  do {
    if (k > max_indel_size) {
      logger.warn("Paramater file '%s' deletion distribution is being "
                  "truncated at max indel length of %d.",
                  param_file_name.c_str(), max_indel_size);
      break;
    }
    _delete_params.increment(k, log(strtod(p,NULL)));
    p = strtok(NULL, "\t");
    k++;
  } while (p);
  
  fix();
}

void MismatchTable::get_indices(const FragHit& f,
                           vector<char>& left_indices,
                           vector<char>& left_seq,
                           vector<char>& left_ref,
                           vector<char>& right_indices,
                           vector<char>& right_seq,
                           vector<char>& right_ref) const {

  const Target& targ = *f.target();
  const Sequence& t_seq_fwd = targ.seq(0);
  const Sequence& t_seq_rev = targ.seq(1);
  
  if (f.left_read()) {
    const ReadHit& read_l = *f.left_read();
    
    left_indices = vector<char>();
    left_seq = vector<char>();
    
    size_t i = 0;  // read index
    size_t j = read_l.left;  // genomic index
    
    vector<Indel>::const_iterator ins = read_l.inserts.begin();
    vector<Indel>::const_iterator del = read_l.deletes.begin();
    
    if (read_l.inserts.size() || read_l.deletes.size()) {
      logger.severe("Indels are not currently supported for eXpress-D.");
    }
    
    size_t cur_seq_bit = 0;
    
    while (i < read_l.seq.length()) {
      if (del != read_l.deletes.end() && del->pos == i) {
        j += del->len;
        del++;
      } else if (ins != read_l.inserts.end() && ins->pos == i) {
        i += ins->len;
        ins++;
      } else {
        size_t cur = read_l.seq[i];
        size_t ref = t_seq_fwd[j];
        if (cur != ref) {
          left_indices.push_back(i);
          if (cur_seq_bit / 8 == left_seq.size()) {
            left_seq.push_back(0);
            left_ref.push_back(0);
          }
          left_seq.back() += cur << (cur_seq_bit % 8);
          left_ref.back() += ref << (cur_seq_bit % 8);
          cur_seq_bit += 2;
        }
        
        i++;
        j++;
      }
    }
  }
  
  if (f.right_read()) {
    const ReadHit& read_r = *f.right_read();
    
    right_indices = vector<char>();
    right_seq = vector<char>();

    size_t r_len = read_r.seq.length();
    size_t i = 0;
    size_t j = targ.length() - read_r.right;
        
    vector<Indel>::const_iterator ins = read_r.inserts.end()-1;
    vector<Indel>::const_iterator del = read_r.deletes.end()-1;

    if (read_r.inserts.size() || read_r.deletes.size()) {
      logger.severe("Indels are not currently supported for eXpress-D.");
    }
    
    size_t cur_seq_bit = 0;
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
        size_t ref = t_seq_rev[j];
        
        if (cur != ref) {
          right_indices.push_back(i);
          if (cur_seq_bit / 8 == right_seq.size()) {
            right_seq.push_back(0);
            right_ref.push_back(0);
          }
          right_seq.back() += cur << (cur_seq_bit % 8);
          right_ref.back() += ref << (cur_seq_bit % 8);
          cur_seq_bit += 2;
        }
        
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

    bool insertion = false;
    bool deletion = false;
    
    vector<Indel>::const_iterator ins = read_l.inserts.begin();
    vector<Indel>::const_iterator del = read_l.deletes.begin();
    
    while (i < read_l.seq.length()) {

      if (del != read_l.deletes.end() && del->pos == i) {
        // Deletion at this position
        ll += _delete_params(del->len);
        j += del->len;
        del++;
        deletion = true;
      } else if (ins != read_l.inserts.end() && ins->pos == i) {
        // Insertion at this position
        ll += _insert_params(ins->len);
        i += ins->len;
        ins++;
        insertion = true;
      } else {
        ll += !insertion * _insert_params(0);
        ll += !deletion * _delete_params(0);
        insertion = false;
        deletion = false;
        
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

    bool insertion = false;
    bool deletion = false;
    
    vector<Indel>::const_iterator ins = read_r.inserts.end()-1;
    vector<Indel>::const_iterator del = read_r.deletes.end()-1;

    while (i < r_len) {
      if (del != read_r.deletes.begin() - 1 && del->pos == r_len - i) {
        ll += _delete_params(del->len);
        j += del->len;
        del--;
        deletion = true;
      } else if (ins != read_r.inserts.begin() - 1
                 && ins->pos + ins->len == r_len - i) {
        ll += _insert_params(ins->len);
        i += ins->len;
        ins--;
        insertion = true;
      } else {
        ll += !insertion * _insert_params(0);
        ll += !deletion * _delete_params(0);
        insertion = false;
        deletion = false;
        
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
          size_t index = prev + ref;
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
  if (mass == LOG_0) {
    return;
  }

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

    bool insertion = false;
    bool deletion = false;
    
    vector<Indel>::const_iterator ins = read_l.inserts.begin();
    vector<Indel>::const_iterator del = read_l.deletes.begin();

    vector<double> joint_probs(NUM_NUCS);

    assert(targ.length() >= f.right());
    while (i < read_l.seq.length()) {
      if (del != read_l.deletes.end() && del->pos == i) {
        _delete_params.increment(del->len, mass);
        j += del->len;
        del++;
        deletion = true;
      } else if (ins != read_l.inserts.end() && ins->pos == i) {
        _insert_params.increment(ins->len, mass);
        i += ins->len;
        ins++;
        insertion = true;
      } else {
        if (!insertion) {
          _insert_params.increment(0, mass);
        }
        if (!deletion) {
          _delete_params.increment(0, mass);
        }
        insertion = false;
        deletion = false;
        
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

    bool insertion = false;
    bool deletion = false;
    
    vector<Indel>::const_iterator ins = read_r.inserts.end() - 1;
    vector<Indel>::const_iterator del = read_r.deletes.end() - 1;

    vector<double> joint_probs(NUM_NUCS);

    while (i < r_len) {
      if (del != read_r.deletes.begin()-1 && del->pos == r_len-i ) {
        _delete_params.increment(del->len, mass);
        j += del->len;
        del--;
        deletion = true;
      } else if (ins != read_r.inserts.begin() - 1 &&
                 ins->pos + ins->len == r_len-i) {
        _insert_params.increment(ins->len, mass);
        i += ins->len;
        ins--;
        insertion = true;
      } else {
        if (!insertion) {
          _delete_params.increment(0, mass);
        }
        if (!deletion) {
          _insert_params.increment(0, mass);
        }
        insertion = false;
        deletion = false;

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
          right_mm[i].increment(index, cur, mass+p);
        }

        i++;
        j++;
      }
      _max_len = max(_max_len, read_r.seq.length());
    }
  }
}

void MismatchTable::fix() {
  for (size_t i = 0; i < max_read_len; i++) {
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
    outfile << k+1 << ":\t";
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
    outfile << k+1 << ":\t";
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
  outfile << ">Insertion Length (0-" << max_indel_size << ")\n";
  for (size_t i = 0; i <= max_indel_size; i++) {
    outfile << scientific << sexp(_insert_params(i))<<"\t";
  }
  outfile<<endl;
  outfile << ">Deletion Length (0-" << max_indel_size << ")\n";
  for (size_t i = 0; i <= max_indel_size; i++) {
    outfile << scientific << sexp(_delete_params(i))<<"\t";
  }
  outfile<<endl;
}

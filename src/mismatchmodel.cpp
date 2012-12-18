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
      _insert_params(1, max_indel_size + 1, 0),
      _delete_params(1, max_indel_size + 1, 0),
      _max_len(0),
      _active(false){
  // Set indel priors
  double indel_p = 1. - pow(EPSILON/2, 1/((double)max_indel_size + 1.));
  double pm = indel_p;
  for(size_t i = 0 ; i <= max_indel_size; ++i) {
    _insert_params.increment(i, log(alpha * pm));
    _delete_params.increment(i, log(alpha * pm));
    pm *= (1 - indel_p);
  }
  assert(approx_eq(sexp(_insert_params.sum(0)), alpha));
  assert(approx_eq(sexp(_delete_params.sum(0)), alpha));
}

MismatchTable::MismatchTable(string param_file_name)
    : _first_read_mm(MAX_READ_LEN, FrequencyMatrix<double>(16, 4, 0)),
      _second_read_mm(MAX_READ_LEN, FrequencyMatrix<double>(16, 4, 0)),
      _insert_params(1, max_indel_size + 1, 0),
      _delete_params(1, max_indel_size + 1, 0),
      _max_len(0),
      _active(true){
  ifstream infile (param_file_name.c_str());
  size_t BUFF_SIZE = 99999;
  char line_buff[BUFF_SIZE];
  if (!infile.is_open()) {
    cerr << "ERROR: Unable to open parameter file '" << param_file_name
         << "'.\n";
    exit(1);
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
    
    if (pos >= MAX_READ_LEN) {
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
  _max_len = min(pos, MAX_READ_LEN);
  if (pos >= MAX_READ_LEN) {
    cerr << "WARNING: First read error distribution of " << pos-1
         << " bases in '" << param_file_name << "' truncated after "
         << MAX_READ_LEN << " bases.\n";
  }
  
  pos = 0;
  infile.getline (line_buff, BUFF_SIZE, '\n');
  
  while ( infile.good() ) {
    infile.getline (line_buff, BUFF_SIZE, '\n');
    
    if (!strncmp(line_buff, ">Insertion Length", 17)) {
      break;
    }
    
    if (pos >= MAX_READ_LEN) {
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
  _max_len = max(_max_len, min(pos, MAX_READ_LEN));
  if (pos >= MAX_READ_LEN) {
    cerr << "WARNING: Second read error distribution of " << pos-1
    << " bases in '" << param_file_name << "' truncated after "
    << MAX_READ_LEN << " bases.\n";
  }
  
  infile.getline (line_buff, BUFF_SIZE, '\n');
  char *p = strtok(line_buff, "\t");
  size_t k = 0;
  do {
    if (k > max_indel_size) {
      cerr << "WARNING: Paramater file '" << param_file_name << "' insertion "
           << "distribution is being truncated at max indel length of "
           << max_indel_size << ".\n";
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
      cerr << "WARNING: Paramater file '" << param_file_name << "' deletion "
           << "distribution is being truncated at max indel length of "
           << max_indel_size << ".\n";
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

    bool insertion, deletion;
    
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

    bool insertion, deletion;
    
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

    bool insertion, deletion;
    
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

    bool insertion, deletion;
    
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

/**
 *  targets.cpp
 *  express
 *
 *  Created by Adam Roberts on 3/20/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 */

#include "main.h"
#include "targets.h"
#include "fld.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include "mapparser.h"
#include "library.h"
#include "tautree.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>

using namespace std;

Target::Target(TargID id, const std::string& name, const std::string& seq,
               bool prob_seq)
   : _id(id),
     _name(name),
     _seq_f(seq, 0, prob_seq),
     _seq_r(_seq_f),
     _ret_params(&_curr_params),
     _uniq_counts(0),
     _tot_counts(0),
     _avg_bias(0),
     _cached_eff_len(log(seq.length())),
     _solvable(false) {
}


void Target::round_reset() {
  _last_params = _curr_params;
  _curr_params = RoundParams();
  _ret_params = &_last_params;
}

double Target::est_effective_length(const FLD* fld, bool with_bias) const {
  double eff_len = LOG_0;

  for(size_t l = fld->min_val(); l <= min(length(), fld->max_val()); l++) {
    eff_len = log_add(eff_len, fld->pmf(l)+log((double)length()-l+1));
  }

  if (with_bias) {
    boost::shared_lock<boost::shared_mutex> lock(_mutex);
    eff_len += _avg_bias;
  }

  return eff_len;
}

double Target::cached_effective_length(bool with_bias) const {
  boost::shared_lock<boost::shared_mutex> lock(_mutex);
  if (with_bias) {
    return _cached_eff_len + _avg_bias;
  }
  return _cached_eff_len;
}

float Target::get_5_bias(size_t pos) const {
  boost::shared_lock<boost::shared_mutex> lock(_mutex);
  if (_5_bias) {
    return _5_bias->at(pos);
  }
  return 0;
}

float Target::get_3_bias(size_t pos) const {
  boost::shared_lock<boost::shared_mutex> lock(_mutex);
  if (_3_bias) {
    return _3_bias->at(pos);
  }
  return 0;
}

void Target::update_target_bias(BiasBoss* bias_table, FLD* fld) {
  boost::unique_lock<boost::shared_mutex> lock(_mutex);

  if (bias_table) {
    if (!_5_bias) {
      _5_bias.reset(new std::vector<float>(length(), 0));
      _3_bias.reset(new std::vector<float>(length(), 0));
    }
    _avg_bias = bias_table->get_target_bias(*_3_bias, *_5_bias, *this);
    assert(!isnan(_avg_bias) && !isinf(_avg_bias));
  }
  _cached_eff_len = est_effective_length(fld, false);
}

TargetTable::TargetTable(const string& targ_fasta_file, bool prob_seqs,
                         const Librarian* libs)
    :  _libs(libs) {
  cout << "Loading target sequences";
  const Library& lib = _libs->curr_lib();
  const TransIndex& targ_index = lib.map_parser->targ_index();
  const TransIndex& targ_lengths = lib.map_parser->targ_lengths();
  if (lib.bias_table) {
        cout << " and measuring bias background";
  }
  cout << "...\n\n";

  size_t num_targs = targ_index.size();
  _targ_map = vector<Target*>(num_targs, NULL);

  boost::unordered_set<string> target_names;


  ifstream infile (targ_fasta_file.c_str());
  string line;
  string seq = "";
  string name = "";
  if (infile.is_open()) {
    while (infile.good()) {
      getline(infile, line, '\n');
      if (line[0] == '>') {
        if (!name.empty()) {
          add_targ(name, seq, prob_seqs, targ_index, targ_lengths);
        }
        name = line.substr(1,line.find(' ')-1);
        if (target_names.count(name)) {
          cerr << "ERROR: Target '" << name << "' is duplicated in the input "
               << "FASTA. Ensure all target names are unique and re-map before "
               << "re-running eXpress\n";
          exit(1);
        }
        target_names.insert(name);
        seq = "";
      } else {
        seq += line;
      }
    }
    if (!name.empty()) {
      add_targ(name, seq, prob_seqs, targ_index, targ_lengths);
    }

    infile.close();
    if (lib.bias_table) {
            lib.bias_table->normalize_expectations();
    }
  } else {
    cerr << "ERROR: Unable to open MultiFASTA file '" << targ_fasta_file
         << "'.\n" ;
    exit(1);
  }

  if (size() == 0) {
    cerr << "ERROR: No targets found in MultiFASTA file '" << targ_fasta_file
         << "'.\n" ;
    exit(1);
  }

  for(TransIndex::const_iterator it = targ_index.begin();
      it != targ_index.end(); ++it) {
    if (!_targ_map[it->second]) {
      cerr << "ERROR: Sequence for target '" << it->first
           << "' not found in MultiFasta file '" << targ_fasta_file << "'.\n";
      exit(1);
    }
  }
}

TargetTable::~TargetTable() {
  foreach( Target* targ, _targ_map) {
    delete targ;
  }
}

vector<double> TargetTable::get_alphas(double alpha,
                                       const AlphaMap* alpha_map) const {
  double alpha_renorm = 1.0;
  if (alpha_map) {
    double alpha_total = 0;
    for (AlphaMap::const_iterator it = alpha_map->begin();
         it != alpha_map->end(); ++it) {
      alpha_total += it->second;
    }
    alpha_renorm = (alpha * alpha_map->size())/alpha_total;
  }
  
  vector<double> target_alphas(_targ_map.size(), LOG_0);
  for (size_t i = 0; i < _targ_map.size(); ++i) {
    if (alpha_map) {
      string name = _targ_map[i]->name();
      
      if (!alpha_map->count(name)) {
        cerr << "ERROR: Target '" << name << "' is was not found in the "
        << "prior parameter file.\n";
        exit(1);
      }

      alpha = alpha_renorm * alpha_map->find(name)->second;
    }
    target_alphas[i] = log(alpha) + _targ_map[i]->cached_effective_length();
  }
  return target_alphas;
}

void TargetTable::add_targ(const string& name, const string& seq, bool prob_seq,
                           const TransIndex& targ_index,
                           const TransIndex& targ_lengths) {
  TransIndex::const_iterator it = targ_index.find(name);
  if (it == targ_index.end()) {
    cerr << "Warning: Target '" << name
         << "' exists in MultiFASTA but not alignment (SAM/BAM) file.\n";
    return;
  }

  if (targ_lengths.find(name)->second != seq.length()) {
    cerr << "ERROR: Target '" << name << "' differs in length between "
         << "MultiFASTA and alignment (SAM/BAM) files ("<< seq.length()
         << " vs. " << targ_lengths.find(name)->second << ").\n";
    exit(1);
  }

  Target* targ = new Target(it->second, name, seq, prob_seq);
  const Library& lib = _libs->curr_lib();
  if (lib.bias_table) {
      (lib.bias_table)->update_expectations(*targ);
  }
  _targ_map[targ->id()] = targ;
  targ->bundle(_bundle_table.create_bundle(targ));
}

Target* TargetTable::get_targ(TargID id) {
    return _targ_map[id];
}

Bundle* TargetTable::merge_bundles(Bundle* b1, Bundle* b2) {
  if (b1 != b2) {
    return _bundle_table.merge(b1, b2);
  }
  return b1;
}

void TargetTable::round_reset() {
  foreach(Target* targ, _targ_map) {
    targ->round_reset();
  }
}

void project_to_polytope(vector<const Target*> bundle_targ,
                         vector<double>& targ_counts, double bundle_counts) {
  vector<bool> polytope_bound(bundle_targ.size(), false);
  while (true) {
    double unbound_counts = 0;
    double bound_counts = 0;
    for (size_t i = 0; i < bundle_targ.size(); ++i) {
      const Target& targ = *bundle_targ[i];
      assert(targ.tot_counts() <= bundle_counts);
      if (targ_counts[i] > targ.tot_counts()) {
        targ_counts[i] = targ.tot_counts();
        polytope_bound[i] = true;
      } else if (targ_counts[i] < targ.uniq_counts()) {
        targ_counts[i] = targ.uniq_counts();
        polytope_bound[i] = true;
      }
      
      if (polytope_bound[i]) {
        bound_counts += targ_counts[i];
      } else {
        unbound_counts += targ_counts[i];
      }
    }
   
    if (approx_eq(unbound_counts + bound_counts, bundle_counts)) {
      return;
    }
	
	  if (unbound_counts == 0) {
      polytope_bound = vector<bool>(bundle_targ.size(), false);
	    unbound_counts = bound_counts;
	    bound_counts = 0;
	  }

    double normalizer = (bundle_counts - bound_counts)/unbound_counts;
    for (size_t i = 0; i < bundle_targ.size(); ++i) {
      if (!polytope_bound[i]) {
        targ_counts[i] *= normalizer;
      }
    }
  }
}

void TargetTable::output_results(string output_dir, size_t tot_counts,
                                 const RangeTauForest* forest, const FLD* fld,
                                 bool output_rdds) {
  FILE * expr_file = fopen((output_dir + "/results.xprs").c_str(), "w");
  ofstream varcov_file;
  ofstream   rdds_file;
  if (output_rdds) {
      rdds_file.open((output_dir + "/  rdds.xprs").c_str());
      rdds_file << "target_id\tposition\tp_value\tref_nuc\tP(A)\tP(C)\tP(G)\t"
               << "P(T)\tobs_A\tobs_C\tobs_G\tobs_T\texp_A\texp_C\texp_G\texp_T"
               << endl;
  }

  fprintf(expr_file, "bundle_id\ttarget_id\tlength\teff_length\ttot_counts\t"
                     "uniq_counts\test_counts\teff_counts\tambig_distr_alpha\t"
                     "ambig_distr_beta\tfpkm\tfpkm_conf_low\tfpkm_conf_high\t"
                     "solvable\n");

  double l_bil = log(1000000000.);
  double l_tot_counts = log((double)tot_counts);

  for (size_t bundle_id = 0; bundle_id < forest->num_children(); ++bundle_id) {
    RangeTauTree& tree = *static_cast<RangeTauTree*>(forest->child(bundle_id));
    double l_bundle_counts = log((double)forest->tree_counts(bundle_id));
    
    if (forest->tree_counts(bundle_id)) {
      vector<double> targ_counts(tree.num_leaves(),0);
      vector<const Target*> bundle_targ(tree.num_leaves(), NULL);
      bool requires_projection = false;

      TauLeafIterator it(&tree);
      size_t i = 0;
      while (*it != NULL) {
        assert(i == (*it)->id() - tree.left());
        bundle_targ[i] = _targ_map[(*it)->id()];
        const Target& targ = *bundle_targ[i];
        targ_counts[i] = sexp(it.tau() + l_bundle_counts);
        assert(!isnan(targ_counts[i]));
        requires_projection |= targ_counts[i] > (double)targ.tot_counts() ||
                               targ_counts[i] < (double)targ.uniq_counts();

        ++i;
        ++it;
      }
                   
      if (tree.num_leaves() > 1 && requires_projection) {
        project_to_polytope(bundle_targ, targ_counts,
                            forest->tree_counts(bundle_id));
      }
       
      // Calculate individual counts and taus
      for (size_t i = 0; i < tree.num_leaves(); ++i) {
        const Target& targ = *bundle_targ[i];
        double l_eff_len = targ.est_effective_length(fld);

        double fpkm_constant = sexp(l_bil - l_eff_len - l_tot_counts);
        double targ_fpkm = targ_counts[i] * fpkm_constant;
        
        double eff_len = sexp(l_eff_len);
        double eff_counts = targ_counts[i] / eff_len * targ.length();
        
        double count_alpha = -1;
        double count_beta = -1;
        double fpkm_lo = -1;
        double fpkm_hi = -1;
        
        fprintf(expr_file, "" SIZE_T_FMT "\t%s\t" SIZE_T_FMT "\t%f\t" SIZE_T_FMT
                          "\t" SIZE_T_FMT "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%c\n",
               bundle_id, targ.name().c_str(), targ.length(), eff_len,
               targ.tot_counts(), targ.uniq_counts(), targ_counts[i],
               eff_counts, count_alpha, count_beta, targ_fpkm, fpkm_lo, fpkm_hi,
               (targ.solvable())?'T':'F');
       
    
         if (output_rdds) {
           const Sequence& targ_seq = targ.seq();
           vector<double> p_vals;
           targ_seq.calc_p_vals(p_vals);
           for (size_t i = 0; i < p_vals.size(); ++i) {
             if (p_vals[i] < 0.05) {
                 rdds_file << targ.name() << "\t" << i << "\t" << p_vals[i]
                          << "\t" << NUCS[targ_seq.get_ref(i)];
               for (size_t nuc=0; nuc < NUM_NUCS; nuc++) {
                   rdds_file << "\t" << sexp(targ_seq.get_prob(i,nuc));
               }
               for (size_t nuc=0; nuc < NUM_NUCS; nuc++) {
                   rdds_file << "\t" << sexp(targ_seq.get_obs(i,nuc));
               }
               for (size_t nuc=0; nuc < NUM_NUCS; nuc++) {
                   rdds_file << "\t" << sexp(targ_seq.get_exp(i,nuc));
               }
                 rdds_file << endl;
             }
           }
         }
      }
    } else {
      for (size_t i = 0; i < tree.num_leaves(); ++i) {
        Target& targ = *_targ_map[i + tree.left()];
        fprintf(expr_file, "" SIZE_T_FMT "\t%s\t" SIZE_T_FMT "\t%f\t%d\t%d\t%f"
                           "\t%f\t%f\t%f\t%f\t%f\t%f\t%c\n",
                bundle_id, targ.name().c_str(), targ.length(),
                sexp(targ.est_effective_length(fld)), 0, 0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 'T');
      }
    }
  }
  fclose(expr_file);
  if (output_rdds) {
      rdds_file.close();
  }
}

void TargetTable::asynch_bias_update(boost::mutex* mutex) {
  BiasBoss* bg_table = NULL;
  boost::scoped_ptr<BiasBoss> bias_table;
  boost::scoped_ptr<FLD> fld;

  bool burned_out_before = false;

  const Library& lib = _libs->curr_lib();

  while(running) {
    if (bg_table) {
       bg_table->normalize_expectations();
    }
    {
      boost::unique_lock<boost::mutex> lock(*mutex);
      if(!fld) {
        fld.reset(new FLD(*(lib.fld)));
      } else {
        *fld = *(lib.fld);
      }
      if (lib.bias_table) {
        BiasBoss& lib_bias_table = *(lib.bias_table);
        if (!bias_table) {
          bias_table.reset(new BiasBoss(lib_bias_table));
        } else {
          lib_bias_table.copy_expectations(*bg_table);
          bg_table->copy_observations(lib_bias_table);
          bias_table.reset(bg_table);
        }
        bg_table = new BiasBoss(0);
      }
      cout << "Synchronized parameter tables.\n";
    }

    if (!edit_detect && burned_out && burned_out_before) {
      break;
    }

    burned_out_before = burned_out;
   
    vector<double> fl_cdf = fld->cmf();

    TauLeafIterator it(lib.tau_forest);
    while (*it != NULL) {
      TauTree* leaf = *it;
      Target* targ = lib.targ_table->get_targ(leaf->id());
      targ->update_target_bias(bias_table.get(), fld.get());
      if (bg_table) {
        assert(false); // These are taus, not rhos!
        bg_table->update_expectations(*targ, it.tau(), fl_cdf);
      }
      ++it;
      if (!running) {
        break;
      }
    }
  }

  if (bg_table) {
     delete bg_table;
  }
}
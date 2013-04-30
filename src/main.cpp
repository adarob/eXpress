/**
 *  main.cpp
 *  express
 *
 *  Created by Adam Roberts on 3/23/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 **/

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "main.h"
#include "bundles.h"
#include "targets.h"
#include "lengthdistribution.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include "mapparser.h"
#include "threadsafety.h"
#include "robertsfilter.h"
#include "directiondetector.h"
#include "library.h"

#ifdef PROTO
  #include "proto/alignments.pb.h"
  #include "proto/targets.pb.h"
  #include <boost/archive/iterators/base64_from_binary.hpp>
  #include <boost/archive/iterators/transform_width.hpp>
  #include <boost/archive/iterators/ostream_iterator.hpp>
#endif

#ifndef WIN32
  #include "update_check.h"
#endif

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

// the forgetting factor parameter controls the growth of the fragment mass
double ff_param = 0.85;

// the burn-in parameter determines how many reads are required before the
// error and bias models are applied to probabilistic assignment
size_t burn_in = 100000;
size_t burn_out = 5000000;
bool burned_out = false;

size_t max_indel_size = 10;

size_t stop_at = 0;

// file location parameters
string output_dir = ".";
string fasta_file_name = "";
string in_map_file_names = "";
string param_file_name = "";
string haplotype_file_name = "";

// intial pseudo-count parameters (non-logged)
double expr_alpha = .005;
double fld_alpha = 1;
double bias_alpha = 1;
double mm_alpha = 1;

size_t bias_model_order = 3;

// fragment length parameters
size_t def_fl_max = 800;
size_t def_fl_mean = 200;
size_t def_fl_stddev = 80;
size_t def_fl_kernel_n = 4;
double def_fl_kernel_p = 0.5;

// option parameters
bool edit_detect = false;
bool error_model = true;
bool bias_correct = true;
bool calc_covar = false;
bool output_align_prob = false;
bool output_align_samp = false;
bool output_running_rounds = false;
bool output_running_reads = false;
size_t num_threads = 2;
size_t num_neighbors = 0;
size_t library_size = 0;

// directional parameters
Direction direction = BOTH;

bool running = true;

// used for multiple rounds of EM
bool first_round = true;
bool last_round = true;
bool batch_mode = false;
bool online_additional = false;
bool both = false;
size_t remaining_rounds = 0;

bool spark_pre = false;

typedef boost::unordered_map<string, double> AlphaMap;
AlphaMap* expr_alpha_map = NULL;

/**
 * Parses an input file of pseudo-count priors for targets.
 * @param in_file path to the input file.
 * @return A pointer to a mapping from target names to prior. Must be deleted.
 */
AlphaMap* parse_priors(string in_file) {
  ifstream ifs(in_file.c_str());
  if (!ifs.is_open()) {
    cerr << "ERROR: Unable to open input priors file '" << in_file << "'.\n" ;
    exit(1);
  }
  AlphaMap* alphas = new AlphaMap();

  string line;

  while(ifs.good()) {
    getline(ifs,line);

    size_t idx = line.find_first_of("\t ");
    if (idx!=string::npos) {
      string name = line.substr(0,idx);
      string val = line.substr(idx+1);
      (*alphas)[name] = atof(val.c_str());
    }
  }
  return alphas;
};

/**
 * Parses argument options and sets variables appropriately.
 * @param ac number of arguments.
 * @param pointer to array of arguments as character arrays.
 * @return True iff there was an error.
 */
bool parse_options(int ac, char ** av) {

  size_t additional_online = 0;
  size_t additional_batch = 0;
  
  po::options_description standard("Standard Options");
  standard.add_options()
  ("help,h", "produce help message")
  ("output-dir,o", po::value<string>(&output_dir)->default_value(output_dir),
   "write all output files to this directory")
#ifdef PROTO
  ("preprocess,D", "run preprocess script for eXpressD")
#endif
  ("frag-len-mean,m", po::value<size_t>(&def_fl_mean)->default_value(def_fl_mean),
   "prior estimate for average fragment length")
  ("frag-len-stddev,s",
   po::value<size_t>(&def_fl_stddev)->default_value(def_fl_stddev),
   "prior estimate for fragment length std deviation")
  ("haplotype-file,H",
   po::value<string>(&haplotype_file_name)->default_value(haplotype_file_name),
   "path to a file containing haplotype pairs")
  ("additional-batch,B",
   po::value<size_t>(&additional_batch)->default_value(additional_batch),
   "number of additional batch EM rounds after initial online round")
  ("additional-online,O",
   po::value<size_t>(&additional_online)->default_value(additional_online),
   "number of additional online EM rounds after initial online round")
  ("output-align-prob",
   "output alignments (sam/bam) with probabilistic assignments")
  ("output-align-samp",
   "output alignments (sam/bam) with sampled assignments")
  ("fr-stranded",
   "accept only forward->reverse alignments (second-stranded protocols)")
  ("rf-stranded",
   "accept only reverse->forward alignments (first-stranded protocols)")
  ("f-stranded",
   "accept only forward single-end alignments (second-stranded protocols)")
  ("r-stranded",
   "accept only reverse single-end alignments (first-stranded protocols)")
  ("no-update-check", "disables automatic check for update via web")
  ;
  
  po::options_description advanced("Advanced Options");
  advanced.add_options()
  ("forget-param,f", po::value<double>(&ff_param)->default_value(ff_param),
   "sets the 'forgetting factor' parameter (0.5 < c <= 1)")
  ("library-size", po::value<size_t>(&library_size),
   "specifies library size for FPKM instead of calculating from alignments")
  ("max-indel-size",
   po::value<size_t>(&max_indel_size)->default_value(max_indel_size),
   "sets the maximum allowed indel size, affecting geometric indel prior")
  ("calc-covar", "calculate and output covariance matrix")
  ("expr-alpha", po::value<double>(&expr_alpha)->default_value(expr_alpha),
   "sets the strength of the prior, per bp")
  ("stop-at", po::value<size_t>(&stop_at)->default_value(stop_at),
   "sets the number of fragments to process, disabled with 0")
  ("burn-out", po::value<size_t>(&burn_out)->default_value(burn_out),
   "sets number of fragments after which to stop updating auxiliary parameters")
  ("no-bias-correct", "disables bias correction")
  ("no-error-model", "disables error modelling")
  ("aux-param-file",
   po::value<string>(&param_file_name)->default_value(param_file_name),
   "path to file containing auxiliary parameters to use instead of learning")
  ;

  string prior_file = "";

  po::options_description hidden("Experimental/Debug Options");
  hidden.add_options()
  ("num-threads,p", po::value<size_t>(&num_threads)->default_value(num_threads),
   "number of threads (>= 2)")
  ("edit-detect","")
  ("single-round", "")
  ("output-running-rounds", "")
  ("output-running-reads", "")
  ("batch-mode","")
  ("both","")
  ("prior-params", po::value<string>(&prior_file)->default_value(""), "")
  ("sam-file", po::value<string>(&in_map_file_names)->default_value(""), "")
  ("fasta-file", po::value<string>(&fasta_file_name)->default_value(""), "")
  ("num-neighbors", po::value<size_t>(&num_neighbors)->default_value(0), "")
  ("bias-model-order",
   po::value<size_t>(&bias_model_order)->default_value(bias_model_order),
   "sets the order of the Markov chain used to model sequence bias")
  ;

  po::positional_options_description positional;
  positional.add("fasta-file",1).add("sam-file",1);

  po::options_description cmdline_options;
  cmdline_options.add(standard).add(advanced).add(hidden);

  bool error = false;
  po::variables_map vm;
  try {
    po::store(po::command_line_parser(ac, av).options(cmdline_options)
              .positional(positional).run(), vm);
  } catch (po::error& e) {
    cerr << "Command-Line Argument Error: "<< e.what() << endl;
    error = true;
  }
  po::notify(vm);

  if (ff_param > 1.0 || ff_param < 0.5) {
    cerr << "Command-Line Argument Error: forget-param/f option must be "
         << "between 0.5 and 1.0\n\n";
    error= true;
  }

  if (fasta_file_name == "") {
    cerr << "Command-Line Argument Error: target sequence fasta file "
         << "required\n\n";
    error = true;
  }

  if (error || vm.count("help")) {
    cerr << "express v" << PACKAGE_VERSION << endl
         << "-----------------------------\n"
         << "File Usage:  express [options] <target_seqs.fa> <hits.(sam/bam)>\n"
         << "Piped Usage: bowtie [options] -S <index> <reads.fq> | express "
         << "[options] <target_seqs.fa>\n\n"
         << "Required arguments:\n"
         << " <target_seqs.fa>     target sequence file in fasta format\n"
         << " <hits.(sam/bam)>     read alignment file in SAM or BAM format\n\n"
         << standard
         << advanced;
    return 1;
  }

  if (param_file_name.size()) {
    burn_in = 0;
    burn_out = 0;
    burned_out = true;
  }
  
  size_t stranded_count = 0;
  if (vm.count("fr-stranded")) {
    direction = FR;
    stranded_count++;
  }
  if (vm.count("rf-stranded")) {
    direction = RF;
    stranded_count++;
  }
  if (vm.count("f-stranded")) {
    direction = F;
    stranded_count++;
  }
  if (vm.count("r-stranded")) {
    direction = R;
    stranded_count++;
  }
  if (stranded_count > 1) {
    cerr << "ERROR: Multiple strandedness flags cannot be "
    << "specified in the same run.\n";
    return 1;
  }
  
  edit_detect = vm.count("edit-detect");
  calc_covar = vm.count("calc-covar");
  bias_correct = !(vm.count("no-bias-correct"));
  error_model = !(vm.count("no-error-model"));
  output_align_prob = vm.count("output-align-prob");
  output_align_samp = vm.count("output-align-samp");
  output_running_rounds = vm.count("output-running-rounds");
  output_running_reads = vm.count("output-running-reads");
  batch_mode = vm.count("batch-mode");
  both = vm.count("both");
  remaining_rounds = max(additional_online, additional_batch);
  spark_pre = vm.count("preprocess");

  if (batch_mode) {
    ff_param = 1;
  }
  
  if (additional_online > 0 && additional_batch > 0) {
    cerr << "ERROR: Cannot add both online and batch rounds.";
    return 1;
  } else if (additional_online > 0) {
    online_additional = true;
  }
  
  if (output_align_prob && output_align_samp) {
    cerr << "ERROR: Cannot output both alignment probabilties and sampled "
         << "alignments.";
    return 1;
  }
  if ((output_align_prob || output_align_samp) && remaining_rounds == 0) {
    cerr << "Warning: It is recommended that at least one additional round "
         << "be used when outputting alignment probabilities or sampled "
         << "alignments. Use the '-B' or '-O' option to enable.";
  }
  
  // We have 1 processing thread and 1 parsing thread always, so we should not
  // count these as additional threads.
  if (num_threads < 2) {
    num_threads = 0;
  }
  num_threads -= 2;
  if (num_threads > 0) {
    num_threads -= edit_detect;
  }
  if (remaining_rounds && in_map_file_names == "") {
    cerr << "ERROR: Cannot process multiple rounds from streaming input.";
    return 1;
  }
  if (remaining_rounds) {
    last_round = false;
  }
  if (prior_file != "") {
    expr_alpha_map = parse_priors(prior_file);
  }

#ifndef WIN32
  if (!vm.count("no-update-check")) {
    check_version(PACKAGE_VERSION);
  }
#endif

  return 0;
}

/**
 * This function writes the current abundance parameters to one file and the
 * auxiliary parameters for each library to a separate file.
 * @param libs a Librarian containing the parameters tables for each Library.
 * @param tot_counts a size_t for the total number of fragments processed thus
          far.
 * @param n an int suffix to add to the output subdirectory. No subdirectory is
 *        used if -1 (default).
 */
void output_results(Librarian& libs, size_t tot_counts, int n=-1) {
  char buff[500];
  string dir = output_dir;
  if (n >= 0) {
    sprintf(buff, "%s/x_%d", output_dir.c_str(), n);
    cerr << "Writing results to " << buff << endl;
    dir = string(buff);
    try {
      fs::create_directories(dir);
    } catch (fs::filesystem_error& e) {
            cerr << e.what() << endl;
            exit(1);
    }
  }
  libs[0].targ_table->output_results(dir, tot_counts, last_round&calc_covar,
                                     last_round&edit_detect);

  for (size_t l = 0; l < libs.size(); l++) {
    if (libs.size() > 1) {
      sprintf(buff, "%s/params.%d.xprs", dir.c_str(), (int)l+1);
    } else {
      sprintf(buff, "%s/params.xprs", dir.c_str());
    }
    ofstream paramfile(buff);
    (libs[l].fld)->append_output(paramfile, "Fragment");
    if (libs[l].mismatch_table) {
      (libs[l].mismatch_table)->append_output(paramfile);
    }
    if (libs[l].bias_table) {
      (libs[l].bias_table)->append_output(paramfile);
    }
    paramfile.close();
  }
}

/**
 * This function handles the probabilistic assignment of multi-mapped reads. The
 * marginal likelihoods are calculated for each mapping, and the mass of the
 * fragment is divided based on the normalized marginals to update the model
 * parameters.
 * @param frag_p pointer to the fragment to probabilistically assign.
 */
void process_fragment(Fragment* frag_p) {
  Fragment& frag = *frag_p;
  const Library& lib = *frag.lib();

  // sort hits to avoid deadlock
  frag.sort_hits();
  double mass_n = frag.mass();

  assert(frag.num_hits());

  vector<double> masses(frag.num_hits(), 0);
  vector<double> variances(frag.num_hits(), 0);
  double total_likelihood = LOG_0;
  double total_mass = LOG_0;
  double total_variance = LOG_0;
  size_t num_solvable = 0;

  // set of locked targets
  boost::unordered_set<const Target*> targ_set;
  boost::unordered_set<const Target*> locked_set;

  // Update bundles and merge in first loop
  Bundle* bundle = frag.hits()[0]->target()->bundle();
  
  // calculate marginal likelihoods and lock targets.
  if (frag.num_hits()>1) {
    for (size_t i = 0; i < frag.num_hits(); ++i) {
      FragHit& m = *frag.hits()[i];
      Target* t = m.target();
      
      bundle = lib.targ_table->merge_bundles(bundle, t->bundle());
      
      if (locked_set.count(t) == 0) {
        t->lock();
        locked_set.insert(t);
      }
      targ_set = locked_set;
      // lock neighbors
      foreach (const Target* neighbor, *m.neighbors()) {
        if (locked_set.count(neighbor) == 0) {
          neighbor->lock();
          locked_set.insert(neighbor);
        }
      }
      m.params()->align_likelihood = t->align_likelihood(m);
      m.params()->full_likelihood = m.params()->align_likelihood +
                                    t->sample_likelihood(first_round,
                                                         m.neighbors());
      masses[i] = t->mass();
      variances[i] = t->mass_var();
      total_likelihood = log_add(total_likelihood, m.params()->full_likelihood);
      total_mass = log_add(total_mass, masses[i]);
      total_variance = log_add(total_variance, variances[i]);
      num_solvable += t->solvable();
    }
  } else {
    FragHit& m = *frag.hits()[0];
    Target* t = m.target();
    t->lock();
    locked_set.insert(t);
    total_likelihood = 0;
    m.params()->align_likelihood = 0;
    m.params()->full_likelihood = 0;
    foreach (const Target* neighbor, *frag.hits()[0]->neighbors()) {
      if (targ_set.count(neighbor) == 0) {
        neighbor->lock();
        locked_set.insert(neighbor);
      }
    }
  }

  if (islzero(total_likelihood)){
    assert(expr_alpha_map);
    cerr << "Warning: Fragment '" << frag.name() << "' has 0 likelihood of "
         << "originating from the transcriptome. Skipping.";
    foreach (const Target* t, locked_set) {
      t->unlock();
    }
    return;
  }

  if (first_round) {
    bundle->incr_counts();
  }
  if (first_round || online_additional) {
    bundle->incr_mass(mass_n);
  }
  
  // normalize marginal likelihoods
  for (size_t i = 0; i < frag.num_hits(); ++i) {
    FragHit& m = *frag[i];
    Target* t  = m.target();
    
    double p = m.params()->full_likelihood-total_likelihood;
    m.params()->posterior = p;
    if (targ_set.size() > 1) {
      double v = log_add(variances[i] - 2*total_mass,
                  total_variance + 2*masses[i] - 4*total_mass);
      t->add_hit(m, v, mass_n);
    } else if (i == 0) {
      t->add_hit(m, LOG_0, mass_n);
    }

    // update parameters
    if (first_round) {
      if (i == 0 || frag[i-1]->target_id() != t->id()) {
        t->incr_counts(targ_set.size() <= 1);
      }
      if (!t->solvable() && num_solvable == frag.num_hits()-1) {
        t->solvable(true);
      }
      if ((!burned_out || edit_detect) && lib.mismatch_table) {
        (lib.mismatch_table)->update(m, p, lib.mass_n);
      }
      if (!burned_out) {
        if (m.pair_status() == PAIRED) {
          (lib.fld)->add_val(m.length(), p+lib.mass_n);
        }
        if (lib.bias_table) {
          (lib.bias_table)->update_observed(m, p+lib.mass_n);
        }
      }
    }
    if (calc_covar && (last_round || online_additional)) {
      double var = 2*mass_n + p + log_sub(LOG_1, p);
      lib.targ_table->update_covar(m.target_id(), m.target_id(), var);
      for (size_t j = i+1; j < frag.num_hits(); ++j) {
        const FragHit& m2 = *frag.hits()[j];
        double p2 = m2.params()->full_likelihood-total_likelihood;
        if (sexp(p2) == 0) {
          continue;
        }
        double covar = 2*mass_n + p + p2;
        lib.targ_table->update_covar(m.target_id(), m2.target_id(), covar);
      }
    }
  }

  foreach (const Target* t, locked_set) {
    t->unlock();
  }
}

/**
 * This function processes Fragments asynchronously. Fragments are popped from
 * a threadsafe input queue, processed, and then pushed onto a threadsafe output
 * queue.
 * @param pts pointer to a struct with the input and output Fragment queues.
 */
void proc_thread(ParseThreadSafety* pts) {
  while (true) {
    Fragment* frag = pts->proc_on.pop();
    if (!frag) {
      break;
    }
    process_fragment(frag);
    pts->proc_out.push(frag);
  }
}

/**
 * This is the driver function for the main processing thread. This function
 * updates the current fragment mass for libraries, dispatches fragments to be
 * processed once they are passed by the parsing thread, outputs intermediate
 * results, and handles additional online rounds.
 * @param libs a struct containing pointers to the parameter tables (bias_table,
 *        mismatch_table, fld) and parser for all libraries being processed.
 * @return The total number of fragments processed.
 */
size_t threaded_calc_abundances(Librarian& libs) {
  cerr << "Processing input fragment alignments...\n";
  boost::scoped_ptr<boost::thread> bias_update;

  size_t n = 1;
  size_t num_frags = 0;
  double mass_n = 0;
  cerr << setiosflags(ios::left);

  // For log-scale output
  size_t i = 1;
  size_t j = 6;

  DirectionDetector dir_detector;
  Fragment* frag;
  
  while (true) {
    // Loop through libraries
    for (size_t l = 0; l < libs.size(); l++) {
      Library& lib = libs[l];
      libs.set_curr(l);
      MapParser& map_parser = *lib.map_parser;
      boost::mutex bu_mut;
      // Used to signal bias update thread
      running = true;
      ParseThreadSafety pts(max((int)num_threads,10));
      boost::thread parse(&MapParser::threaded_parse, &map_parser, &pts,
                          stop_at, num_neighbors);
      vector<boost::thread*> thread_pool;
      RobertsFilter frags_seen;

      burned_out = lib.n >= burn_out;
      while(true) {
        if (lib.n == burn_in) {
          bias_update.reset(new boost::thread(&TargetTable::asynch_bias_update,
                                              lib.targ_table, &bu_mut));
          if (lib.mismatch_table) {
            (lib.mismatch_table)->activate();
          }
        }
        if (lib.n == burn_out) {
          (lib.mismatch_table)->fix();
          burned_out = true;
        }
        // Start threads once aux parameters are burned out
        if (burned_out && num_threads && thread_pool.size() == 0) {
          thread_pool = vector<boost::thread*>(num_threads);
          for (size_t k = 0; k < thread_pool.size(); k++) {
            thread_pool[k] = new boost::thread(proc_thread, &pts);
          }
        }

        // Pop next parsed fragment and set mass
        frag = pts.proc_in.pop();
        if (frag) {
          frag->mass(mass_n);
          dir_detector.add_fragment(frag);
        }

        // Test that we have not already seen this fragment
        if (frag && first_round && frags_seen.test_and_push(frag->name())) {
          cerr << "ERROR: Alignments are not properly sorted. Read '"
              << frag->name()
              << "' has alignments which are non-consecutive.\n";
          exit(1);
        }

        // If multi-threaded and burned out, push to the processing queue
        if (num_threads && burned_out) {
          // If no more fragments, send stop signal (NULL) to processing threads
          if (!frag) {
            for (size_t k = 0; k < thread_pool.size(); ++k) {
              pts.proc_on.push(NULL);
            }
            break;
          }
          pts.proc_on.push(frag);
        } else {
          if (!frag) {
            break;
          }
          {
            // Block the bias update thread from updating the paramater tables
            // during processing. We don't need to do this during multi-threaded
            // processing since the parameters are burned out before we start
            // the threads.
            boost::unique_lock<boost::mutex> lock(bu_mut);
            process_fragment(frag);
            pts.proc_out.push(frag);
          }
        }

        // Output intermediate results, if necessary
        if (output_running_reads && n == i*pow(10.,(double)j)) {
          boost::unique_lock<boost::mutex> lock(bu_mut);
          output_results(libs, n, (int)n);
          if (i++ == 9) {
            i = 1;
            j++;
          }
        }
        num_frags++;

        // Output progress
        if (num_frags % 1000000 == 0) {
          cerr << "Fragments Processed (" << lib.in_file_name << "): "
               << setw(9) << num_frags << "\t Number of Bundles: "
               << lib.targ_table->num_bundles() << endl;
          dir_detector.report_if_improper_direction();
        }

        n++;
        lib.n++;
        mass_n += ff_param*log((double)n-1) - log(pow(n,ff_param) - 1);
        lib.mass_n += ff_param*log((double)lib.n-1) -
                      log(pow(lib.n,ff_param) - 1);
      }

      // Signal bias update thread to stop
      running = false;

      parse.join();
      foreach(boost::thread* t, thread_pool) {
        t->join();
      }

      if (bias_update) {
        bias_update->join();
        bias_update.reset(NULL);
      }
    }

    if (online_additional && remaining_rounds--) {
      if (output_running_rounds) {
        output_results(libs, n, (int)remaining_rounds);
      }

      cerr << remaining_rounds << " remaining rounds." << endl;
      first_round = false;
      last_round = (remaining_rounds==0 && !both);
      for (size_t l = 0; l < libs.size(); l++) {
        libs[l].map_parser->write_active(last_round);
        libs[l].map_parser->reset_reader();
      }
      num_frags = 0;
    } else {
      break;
    }
  }

  cerr << "COMPLETED: Processed " << num_frags
       << " mapped fragments, targets are in "
       << libs[0].targ_table->num_bundles() << " bundles\n";

  return num_frags;
}

/**
 * The main function instantiates the library parameter tables and parsers,
 * calls the processing function, and outputs the results. Also handles
 * additional batch rounds.
 */
int estimation_main() {
  
  if (output_dir != ".") {
    try {
      fs::create_directories(output_dir);
    } catch (fs::filesystem_error& e) {
      cerr << e.what() << endl;
    }
  }
  
  if (!fs::exists(output_dir)) {
    cerr << "ERROR: cannot create directory " << output_dir << ".\n";
    exit(1);
  }
  
  // Parse input file names and instantiate Libray structs.
  vector<string> file_names;
  char buff[999];
  strcpy(buff, in_map_file_names.c_str());
  char * pch = strtok (buff,",");
  while (pch != NULL) {
    file_names.push_back(pch);
    pch = strtok (NULL, ",");
  }
  if (file_names.size() == 0) {
    file_names.push_back("");
  }
  Librarian libs(file_names.size());
  for (size_t i = 0; i < file_names.size(); ++i) {
    char out_map_file_name[500] = "";
    if (output_align_prob) {
      sprintf(out_map_file_name, "%s/hits.%d.prob",
              output_dir.c_str(), (int)i+1);
    }
    if (output_align_samp) {
      sprintf(out_map_file_name, "%s/hits.%d.samp",
              output_dir.c_str(), (int)i+1);
    }
    
    libs[i].in_file_name = file_names[i];
    libs[i].out_file_name = out_map_file_name;
    libs[i].map_parser.reset(new MapParser(&libs[i], last_round));

    if (param_file_name.size()) {
      libs[i].fld.reset(new LengthDistribution(param_file_name, "Fragment"));
      libs[i].mismatch_table.reset((error_model) ?
                                              new MismatchTable(param_file_name)
                                              : NULL);
      libs[i].bias_table.reset((bias_correct) ? new BiasBoss(bias_model_order,
                                                         param_file_name):NULL);
    } else {
      libs[i].fld.reset(new LengthDistribution(fld_alpha, def_fl_max,
                                               def_fl_mean, def_fl_stddev,
                                               def_fl_kernel_n,
                                               def_fl_kernel_p));
      libs[i].mismatch_table.reset((error_model) ? new MismatchTable(mm_alpha)
                                                   :NULL);
      libs[i].bias_table.reset((bias_correct) ? new BiasBoss(bias_model_order,
                                                    bias_alpha):NULL);
    }
    if (i > 0 &&
        (libs[i].map_parser->targ_index() != libs[i-1].map_parser->targ_index()
         || libs[i].map_parser->targ_lengths() !=
         libs[i-1].map_parser->targ_lengths())) {
          cerr << "ERROR: Alignment file headers do not match for '"
          << file_names[i-1] << "' and '" << file_names[i] << "'.";
          exit(1);
        }
  }
  
  boost::shared_ptr<TargetTable> targ_table(
                                  new TargetTable(fasta_file_name,
                                                  haplotype_file_name,
                                                  edit_detect,
                                                  param_file_name.size(),
                                                  expr_alpha, expr_alpha_map,
                                                  &libs));
  size_t max_target_length = 0;
  for(size_t tid=0; tid < targ_table->size(); tid++) {
    max_target_length = max(max_target_length,
                            targ_table->get_targ(tid)->length());
  }

  for (size_t i = 0; i < libs.size(); ++i) {
    libs[i].targ_table = targ_table;
    if (bias_correct) {
      libs[i].bias_table->copy_expectations(*(libs.curr_lib().bias_table));
    }
  }
  double num_targ = (double)targ_table->size();
  
  if (calc_covar && (double)SSIZE_MAX < num_targ*(num_targ+1)) {
    cerr << "Warning: Your system is unable to represent large enough values "
         << "for efficiently hashing target pairs. Covariance calculation "
         << "will be disabled.\n";
    calc_covar = false;
  }
  
  if (batch_mode) {
    targ_table->round_reset();
  }
  
  size_t tot_counts = threaded_calc_abundances(libs);
  if (library_size) {
    tot_counts = library_size;
  }
  
  if (!burned_out && bias_correct && param_file_name == "") {
    cerr << "Warning: Not enough fragments observed to accurately learn bias "
         << "paramaters. Either disable bias correction (--no-bias-correct) or "
         << "provide a file containing auxiliary parameters "
         << "(--aux-param-file).\n";
  }
  
  if (both) {
    remaining_rounds = 1;
    online_additional = false;
  }
  
  if (remaining_rounds) {
    targ_table->masses_to_counts();
  }
  
  targ_table->round_reset();
  ff_param = 1.0;

  first_round = false;
  
  while (!last_round) {
    if (output_running_rounds) {
      output_results(libs, tot_counts, (int)remaining_rounds);
    }
    remaining_rounds--;
    cerr << "\nRe-estimating counts with additional round of EM ("
         << remaining_rounds << " remaining)...\n";
    last_round = (remaining_rounds == 0);
    for (size_t l = 0; l < libs.size(); l++) {
      libs[l].map_parser->write_active(last_round);
      libs[l].map_parser->reset_reader();
    }
    tot_counts = threaded_calc_abundances(libs);
    if (library_size) {
      tot_counts = library_size;
    }
    targ_table->round_reset();
  }
  
	cerr << "Writing results to file...\n";
  output_results(libs, tot_counts);
  cerr << "Done\n";
  
  return 0;
}

#ifdef PROTO
inline string base64_encode(const string& to_encode) {
  using namespace boost::archive::iterators;
  typedef base64_from_binary<transform_width<string::const_iterator,6,8> > it_base64_t;
  unsigned int writePaddChars = (3-to_encode.length()%3)%3;
  string base64(it_base64_t(to_encode.begin()), it_base64_t(to_encode.end()));
  base64.append(writePaddChars,'=');
  return base64;
}

int preprocess_main() {
  try {
    fs::create_directories(output_dir);
  } catch (fs::filesystem_error& e) {
      cerr << e.what() << endl;
  }
  
  if (!fs::exists(output_dir)) {
    cerr << "ERROR: cannot create directory " << output_dir << ".\n";
    exit(1);
  }
  
  Librarian libs(1);
  Library& lib = libs[0];
  lib.in_file_name = in_map_file_names;
  lib.out_file_name = "";
  
  lib.map_parser.reset(new MapParser (&lib, false));
  lib.fld.reset(new LengthDistribution(0, 0, 0, 1, 2, 0));
  MarkovModel bias_model(3, 21, 21, 0);
  MismatchTable mismatch_table(0);
  lib.targ_table.reset(new TargetTable (fasta_file_name, "", 0, 0, NULL, NULL,
                                        &libs));
  
  cerr << "Converting targets to Protocol Buffers...\n";
  fstream targ_out((output_dir + "/targets.pb").c_str(),
                   ios::out | ios::trunc);
  string out_buff;
  proto::Target target_proto;
  for (TargID id = 0; id < lib.targ_table->size(); ++id) {
    target_proto.Clear();
    Target& targ = *lib.targ_table->get_targ(id);
    target_proto.set_name(targ.name());
    target_proto.set_id((unsigned int)targ.id());
    target_proto.set_length((unsigned int)targ.length());
    target_proto.set_seq(targ.seq(0).serialize());
    
    target_proto.SerializeToString(&out_buff);
    targ_out << base64_encode(out_buff) << endl;
  }
  targ_out.close();
  
  cerr << "Converting fragment alignments to Protocol Buffers..\n";
  ostream frag_out(cout.rdbuf());
  
  size_t num_frags = 0;
  cerr << setiosflags(ios::left);
  Fragment* frag;
  
  ParseThreadSafety pts(10);
  boost::thread parse(&MapParser::threaded_parse, lib.map_parser.get(), &pts,
                      stop_at, 0);
  RobertsFilter frags_seen;
  proto::Fragment frag_proto;
  while(true) {
    frag_proto.Clear();
    
    // Pop next parsed fragment and set mass
    frag = pts.proc_in.pop();
    
    if (!frag) {
      break;
    }
    
    // Test that we have not already seen this fragment
    if (frags_seen.test_and_push(frag->name())) {
      cerr << "ERROR: Alignments are not properly sorted. Read '"
      << frag->name() << "' has alignments which are non-consecutive.\n";
      exit(1);
    }
    
    frag_proto.set_paired(frag->paired());
    for (size_t i = 0; i < frag->num_hits(); ++i) {
      FragHit& fh = *(*frag)[i];
      proto::FragmentAlignment& align_proto = *frag_proto.add_alignments();
      align_proto.set_target_id((unsigned int)fh.target_id());
      align_proto.set_length((unsigned int)fh.length());
      
      vector<char> left_mm_indices;
      vector<char> left_mm_seq;
      vector<char> right_mm_indices;
      vector<char> right_mm_seq;
      
      mismatch_table.get_indices(fh, left_mm_indices, left_mm_seq,
                                 right_mm_indices, right_mm_seq);
      
      ReadHit* read_l = fh.left_read();
      if (read_l) {
        proto::ReadAlignment& read_proto = *align_proto.mutable_read_l();
        read_proto.set_first(read_l->first);
        read_proto.set_left_pos(read_l->left);
        read_proto.set_right_pos(read_l->right-1);
        read_proto.set_mismatch_indices(string(left_mm_indices.begin(),
                                               left_mm_indices.end()));
        read_proto.set_mismatch_nucs(string(left_mm_seq.begin(),
                                            left_mm_seq.end()));
      }
      
      ReadHit* read_r = fh.right_read();
      if (read_r) {
        proto::ReadAlignment& read_proto = *align_proto.mutable_read_r();
        read_proto.set_first(read_r->first);
        read_proto.set_left_pos(read_r->left);
        read_proto.set_right_pos(read_r->right-1);
        read_proto.set_mismatch_indices(string(right_mm_indices.begin(),
                                               right_mm_indices.end()));
        read_proto.set_mismatch_nucs(string(right_mm_seq.begin(),
                                            right_mm_seq.end()));
      }
    }
    frag_proto.SerializeToString(&out_buff);
    frag_out << base64_encode(out_buff) << endl;
    
    pts.proc_out.push(frag);
    
    num_frags++;
    
    // Output progress
    if (num_frags % 1000000 == 0) {
      cerr << "Fragments Processed: " << setw(9) << num_frags << endl;
    }
  }
  
  parse.join();
  
  return 0;
}
#endif


int main (int argc, char ** argv)
{

  srand((unsigned int)time(NULL));
  int parse_ret = parse_options(argc,argv);
  if (parse_ret) {
    return parse_ret;
  }
  
#ifdef PROTO
  if (spark_pre) {
    return preprocess_main();
  }
#endif
  
  return estimation_main();
}

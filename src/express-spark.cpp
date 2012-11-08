//
//  main.cpp
//  express-spark
//
//  Created by Adam Roberts on 10/15/12.
//
//

#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "targets.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include "mapparser.h"
#include "threadsafety.h"
#include "robertsfilter.h"
#include "library.h"
#include "main.h"

#ifndef WIN32
#include "update_check.h"
#endif

#include "proto/alignments.pb.h"
#include "proto/targets.pb.h"

#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/archive/iterators/ostream_iterator.hpp>

using namespace std;

namespace po = boost::program_options;

string target_output = "targets.pb";
string in_map_filename = "";
string fasta_filename = "";
size_t stop_at = 0;

bool edit_detect = false;
bool running = true;
bool burned_out = true;
size_t max_indel_size = 10;
Direction direction = BOTH;

/**
 * Parses argument options and sets variables appropriately.
 * @param ac number of arguments.
 * @param pointer to array of arguments as character arrays.
 * @return True iff there was an error.
 */
bool parse_options(int ac, char ** av) {
  po::options_description standard("Standard Options");
  standard.add_options()
  ("help,h", "produce help message")
  ("target-output,o",
   po::value<string>(&target_output)->default_value(target_output),
   "file to write target protobufs to")
  ("fr-stranded",
   "accept only forward->reverse alignments (second-stranded protocols)")
  ("rf-stranded",
   "accept only reverse->forward alignments (first-stranded protocols)")
  ;
  
  po::options_description advanced("Advanced Options");
  advanced.add_options()
  ("max-indel-size",
   po::value<size_t>(&max_indel_size)->default_value(max_indel_size),
   "sets the maximum allowed indel size, affecting geometric indel prior")
  ("stop-at", po::value<size_t>(&stop_at)->default_value(stop_at),
   "sets the number of fragments to process, disabled with 0")
  ;
  
  po::options_description hidden("Hidden Options");
  hidden.add_options()
  ("sam-file", po::value<string>(&in_map_filename)->default_value(""), "")
  ("fasta-file", po::value<string>(&fasta_filename)->default_value(""), "")
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
  
  if (fasta_filename == "") {
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
  
  if (vm.count("fr-stranded")) {
    direction = FR;
  }
  
  if (vm.count("rf-stranded")) {
    if (direction != BOTH) {
      cerr << "ERROR fr-stranded and rf-stranded flags cannot both be "
      << "specified in the same run.\n";
      return 1;
    }
    direction = RF;
  }
  
#ifndef WIN32
  if (!vm.count("no-update-check")) {
    check_version(PACKAGE_VERSION);
  }
#endif
  
  return 0;
}

inline string base64_encode(const string& to_encode) {
  using namespace boost::archive::iterators;
  typedef base64_from_binary<transform_width<string::const_iterator,6,8> > it_base64_t;
  unsigned int writePaddChars = (3-to_encode.length()%3)%3;
  string base64(it_base64_t(to_encode.begin()), it_base64_t(to_encode.end()));
  base64.append(writePaddChars,'=');
  return base64;
}

int preprocess_main() {
  Librarian libs(1);
  Library& lib = libs[0];
  lib.in_file_name = in_map_filename;
  lib.out_file_name = "";
  lib.bias_table = NULL;
  MapParser map_parser(&lib, false);
  lib.map_parser = &map_parser;
  lib.fld = new FLD(0, 0, 0, 1);
  MarkovModel bias_model(3, 21, 21, 0);
  MismatchTable mismatch_table(0);
  TargetTable targ_table(fasta_filename, 0, 0,
                         NULL, &libs);
  lib.targ_table = &targ_table;
  
  cerr << "Converting targets to Protocol Buffers...\n";
  fstream targ_out(target_output.c_str(),
                   ios::out | ios::trunc);
  string out_buff;
  proto::Target target_proto;
  for (TargID id = 0; id < targ_table.size(); ++id) {
    target_proto.Clear();
    Target& targ = *targ_table.get_targ(id);
    target_proto.set_name(targ.name());
    target_proto.set_id((unsigned int)targ.id());
    target_proto.set_length((unsigned int)targ.length());
    
    vector<char> bias_indices_l = bias_model.get_indices(targ.seq(0));
    target_proto.set_bias_indices_l(string(bias_indices_l.begin(),
                                           bias_indices_l.end()));
    
    vector<char> bias_indices_r = bias_model.get_indices(targ.seq(1));
    target_proto.set_bias_indices_r(string(bias_indices_r.begin(),
                                           bias_indices_r.end()));
    
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
  boost::thread parse(&MapParser::threaded_parse, &map_parser, &pts,
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
    string first_l_serial;
    string first_r_serial;
    for (size_t i = 0; i < frag->num_hits(); ++i) {
      FragHit& fh = *(*frag)[i];
      proto::FragmentAlignment& align_proto = *frag_proto.add_alignments();
      align_proto.set_target_id((unsigned int)fh.target_id());
      align_proto.set_length((unsigned int)fh.length());
      
      vector<char> left_mm_indices;
      vector<char> right_mm_indices;
      
      mismatch_table.get_indices(fh, left_mm_indices, right_mm_indices);
      
      ReadHit* read_l = fh.left_read();
      if (read_l) {
        proto::ReadAlignment& read_proto = *align_proto.mutable_read_l();
        read_proto.set_first(read_l->first);
        read_proto.set_error_indices(string(left_mm_indices.begin(),
                                            left_mm_indices.end()));
        vector<char> bias_indices;
        size_t start_pos = bias_model.get_indices(fh.target()->seq(0),
                                                  (int)read_l->left,
                                                  bias_indices);
        read_proto.set_bias_start_pos((unsigned int)start_pos);
        read_proto.set_bias_indices(string(bias_indices.begin(),
                                           bias_indices.end()));
        if (i == 0) {
          read_proto.SerializeToString(&first_l_serial);
        } else {
          string serial;
          read_proto.SerializeToString(&serial);
          if (serial == first_l_serial) {
            read_proto.Clear();
          }
        }
      }
      ReadHit* read_r = fh.right_read();
      if (read_r) {
        proto::ReadAlignment& read_proto = *align_proto.mutable_read_r();
        read_proto.set_first(read_r->first);
        read_proto.set_error_indices(string(right_mm_indices.begin(),
                                            right_mm_indices.end()));
        vector<char> bias_indices;
        size_t start_pos = bias_model.get_indices(fh.target()->seq(1),
                                                  (int)(fh.target()->length()
                                                        - read_r->right),
                                                  bias_indices);
        read_proto.set_bias_start_pos((unsigned int)start_pos);
        read_proto.set_bias_indices(string(bias_indices.begin(),
                                           bias_indices.end()));
        
        if (i == 0) {
          read_proto.SerializeToString(&first_r_serial);
        } else {
          string serial;
          read_proto.SerializeToString(&serial);
          if (serial == first_r_serial) {
            read_proto.Clear();
          }
        }
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


int main(int argc, char * argv[])
{
  int parse_ret = parse_options(argc, argv);
  if (parse_ret) {
    return parse_ret;
  }
  
  return preprocess_main();
}


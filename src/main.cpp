//
//  main.cpp
//  expressionline2
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#define PACKAGE_VERSION "INTERNAL"

#include <boost/filesystem.hpp>
#include <getopt.h>
#include <boost/thread.hpp>
#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>
#include "transcripts.h"
#include "fld.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"

using namespace std;

string output_dir = ".";
double expr_alpha = .001;
double fld_alpha = 1;
double bias_alpha = 1;
double mm_alpha = 1;

int def_fl_max = 800;
int def_fl_mean = 200;
int def_fl_stddev = 80;
bool user_provided_fld = false;

bool bias_correct = true;
bool upper_quart_norm = false;

void print_usage()
{
    
    fprintf(stderr, "cuffexpress v%s\n", PACKAGE_VERSION);
//    fprintf(stderr, "linked against Boost version %d\n", BOOST_VERSION);
    fprintf(stderr, "-----------------------------\n"); 
    fprintf(stderr, "File Usage:   cuffexpress [options] <transcripts.fasta> <hits.sam>\n");
    fprintf(stderr, "Piped Usage:  bowtie [options] -S <reads.fastq> | cuffexpress [options] <transcripts.fasta>\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, " -o/--output-dir       write all output files to this directory              [ default:     ./ ]\n");
    fprintf(stderr, " -m/--frag-len-mean    average fragment length (unpaired reads only)         [ default:    200 ]\n");
    fprintf(stderr, " -s/--frag-len-std-dev fragment length std deviation (unpaired reads only)   [ default:     80 ]\n");
    fprintf(stderr, " --upper-quartile-norm use upper-quartile normalization                      [ default:  FALSE ]\n");
    fprintf(stderr, " --no-bias-correct     do not correct for fragment bias                      [ default:  FALSE ]\n");
}

#define OPT_UPPER_QUARTILE_NORM 200
#define OPT_NO_BIAS_CORRECT 201

const char *short_options = "o:m:s";

static struct option long_options[] = {
    // general options
    {"output-dir",              required_argument,		 0,			 'o'},
    {"frag-len-mean",			required_argument,       0,          'm'},
    {"frag-len-std-dev",		required_argument,       0,          's'},
    {"upper-quartile-norm",     no_argument,             0,          OPT_UPPER_QUARTILE_NORM},
    {"no-bias-correct",         no_argument,             0,          OPT_NO_BIAS_CORRECT},
    {0, 0, 0, 0} // terminator
};

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
int parseInt(int lower, const char *errmsg, void (*print_usage)()) {
    long l;
    char *endPtr= NULL;
    l = strtol(optarg, &endPtr, 10);
    if (endPtr != NULL) {
        if (l < lower) {
            cerr << errmsg << endl;
            print_usage();
            exit(1);
        }
        return (int32_t)l;
    }
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}

int parse_options(int argc, char** argv)
{
    int option_index = 0;
    int next_option;
	
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
        switch (next_option) {
			case -1:     /* Done with options. */
				break;
			case 'm':
                user_provided_fld = true;
				def_fl_mean = (uint32_t)parseInt(0, "-m/--frag-len-mean arg must be at least 0", print_usage);
				break;
			case 's':
                user_provided_fld = true;
                def_fl_stddev = (uint32_t)parseInt(0, "-s/--frag-len-std-dev arg must be at least 0", print_usage);
				break;
            case 'o':
                output_dir = optarg;
				break;
			case OPT_UPPER_QUARTILE_NORM:
                upper_quart_norm = true;
                break;
            case OPT_NO_BIAS_CORRECT:
                bias_correct = false;
                break;
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
    return 0;
}

double log_sum(double x, double y)
{
    return x+log(1+exp(y-x));
}

void process_fragment(Fragment* frag_p, TranscriptTable* trans_table, FLD* fld, BiasBoss* bias_table, MismatchTable* mismatch_table)
{
    Fragment& frag = *frag_p;
    
    if (frag.num_maps()==0)
    {
        delete frag_p;
        return;
    }
    if (frag.num_maps()==1)
    {
        const FragMap& m = *frag.maps()[0];
        Transcript* t  = trans_table->get_trans(m.trans_id);
        if (!t)
        {
            cerr << "ERROR: Transcript " << m.name << " not found in reference fasta.";
            exit(1);
        }
        double mass = 1.0;
        t->add_mass(mass);
        fld->add_val(m.length(), mass);
        if (bias_table)
            bias_table->update_observed(m, *t, mass);
        mismatch_table->update(m, *t, mass);
        delete frag_p;
        return;
    }
    
    vector<double> likelihoods(frag.num_maps());
    vector<Transcript*> transcripts(frag.num_maps());
    double total_likelihood = 0.0;
    for(size_t i = 0; i < frag.num_maps(); ++i)
    {
        const FragMap& m = *frag.maps()[i];
        Transcript* t = trans_table->get_trans(m.trans_id);
        transcripts[i] = t;
        assert(t->id() == m.trans_id);
        //t->update_transcript_bias();
        likelihoods[i] = t->log_likelihood(m);
        if (!likelihoods[i])
            continue;
        total_likelihood = (i) ? log_sum(total_likelihood, likelihoods[i]):likelihoods[i];
    }
    if (total_likelihood == 0)
    {
        delete frag_p;
        return;
    }
    double total_mass = 0.0;
    for(size_t i = 0; i < frag.num_maps(); ++i)
    {
        if (!likelihoods[i])
            continue;
        const FragMap& m = *frag.maps()[i];
        Transcript* t  = transcripts[i];
        double mass = exp(likelihoods[i]-total_likelihood);
        assert(!isnan(mass));
        t->add_mass(mass);
        if (bias_table)
            bias_table->update_observed(m, *t, mass);
        mismatch_table->update(m, *t, mass);
        fld->add_val(m.length(), mass);
        total_mass += mass;
    }
    assert(abs(total_mass-1.0) < .0001);
    delete frag_p;
}

void threaded_calc_abundances(MapParser& map_parser, TranscriptTable* trans_table, FLD* fld, BiasBoss* bias_table, MismatchTable* mismatch_table)
{
    ParseThreadSafety ts;
    ts.proc_lk.lock();
    ts.parse_lk.lock();
    boost::thread parse(&MapParser::threaded_parse, &map_parser, &ts);
    
    size_t count = 0;
    ofstream running_expr_file((output_dir + "/running.expr").c_str());
    trans_table->output_header(running_expr_file);
    while(true)
    {
        ts.proc_lk.lock();
        Fragment* frag = ts.next_frag;
        ts.parse_lk.unlock();
        if (!frag)
        {
            ts.proc_lk.unlock();
            return;
        }
//        count += frag->num_maps();
//        if (count % 100000 < frag->num_maps() )
//            cout<< count << "\n";

        count += 1;
        if (count % 1 == 0)
        {
	  //    cout<< count << "\n";
          running_expr_file << count << '\t';  
	  trans_table->output_current(running_expr_file);
        }
        process_fragment(frag, trans_table, fld, bias_table, mismatch_table);
    }
    running_expr_file.close();
}

void calc_abundances(MapParser& map_parser, TranscriptTable* trans_table, FLD* fld, BiasBoss* bias_table, MismatchTable* mismatch_table)
{
    bool fragments_remain = true;
    size_t count = 0;
    while(fragments_remain)
    {
        Fragment* frag = new Fragment();
        fragments_remain = map_parser.next_fragment(*frag);
        count += frag->num_maps();
        if (count % 100000 < frag->num_maps() )
            cout<< count << "\n";
        process_fragment(frag, trans_table, fld, bias_table, mismatch_table);
    }
}

int main (int argc, char ** argv)
{
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
	
    if(optind >= argc)
    {
        print_usage();
        return 1;
    }
    
    if (output_dir != ".")
    {
        boost::filesystem::create_directories(output_dir);
    }
    if (!boost::filesystem::exists(output_dir))
    {
        cerr << "Error: cannot create directory " << output_dir << ".\n";
        exit(1);
    }
    
    string trans_fasta_file_name = argv[optind++];
    string sam_hits_file_name = (optind >= argc) ? "" : argv[optind++];
        
    FLD fld(fld_alpha, def_fl_max);
    BiasBoss* bias_table = (bias_correct) ? new BiasBoss(bias_alpha):NULL;
    MismatchTable mismatch_table(mm_alpha);
    MapParser map_parser(sam_hits_file_name);
    TranscriptTable trans_table(trans_fasta_file_name, expr_alpha, &fld, bias_table, &mismatch_table);

    if (bias_table)
        boost::thread bias_update(&TranscriptTable::threaded_bias_update, &trans_table);
    threaded_calc_abundances(map_parser, &trans_table, &fld, bias_table, &mismatch_table);
    
    fld.output(output_dir);
    mismatch_table.output(output_dir);
    trans_table.output_expression(output_dir);
    return 0;
}

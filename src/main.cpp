//
//  main.cpp
//  express
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#define PACKAGE_VERSION "INTERNAL"

#include <boost/math/distributions/geometric.hpp>
#include <boost/filesystem.hpp>
#include <getopt.h>
#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include "main.h"
#include "transcripts.h"
#include "fld.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include "mapparser.h"

using namespace std;

// the forgetting factor parameter controls the growth of the fragment mass
double ff_param = 1;

string output_dir = ".";

// intial pseudo-count parameters (non-logged)
double expr_alpha = .001;
double fld_alpha = 1;
double bias_alpha = 1;
double mm_alpha = 1;

// fragment length parameters
int def_fl_max = 800;
int def_fl_mean = 200;
int def_fl_stddev = 60;
bool user_provided_fld = false;

// option parameters
bool bias_correct = true;
bool upper_quart_norm = false;
bool calc_covar = false;
bool vis = false;
bool in_between = false;

bool running = true;

void print_usage()
{
    
    fprintf(stderr, "express v%s\n", PACKAGE_VERSION);
//    fprintf(stderr, "linked against Boost version %d\n", BOOST_VERSION);
    fprintf(stderr, "-----------------------------\n"); 
    fprintf(stderr, "File Usage:   cuffexpress [options] <transcripts.fasta> <hits.sam>\n");
    fprintf(stderr, "Piped Usage:  bowtie [options] -S <reads.fastq> | cuffexpress [options] <transcripts.fasta>\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, " -o/--output-dir       write all output files to this directory              [ default:     ./ ]\n");
    fprintf(stderr, " -f/--forget-param     forgetting factor exponent parameter                  [ default:    1.0 ]\n");
    fprintf(stderr, " -m/--frag-len-mean    average fragment length (unpaired reads only)         [ default:    200 ]\n");
    fprintf(stderr, " -s/--frag-len-std-dev fragment length std deviation (unpaired reads only)   [ default:     80 ]\n");
    fprintf(stderr, " --upper-quartile-norm use upper-quartile normalization                      [ default:  FALSE ]\n");
    fprintf(stderr, " --calc-covar          calculate covariance matrix                           [ default:  FALSE ]\n");
    fprintf(stderr, " --no-bias-correct     do not correct for fragment bias                      [ default:  FALSE ]\n");
}

#define OPT_UPPER_QUARTILE_NORM 200
#define OPT_NO_BIAS_CORRECT 201
#define OPT_IN_BETWEEN 202
#define OPT_CALC_COVAR 203

const char *short_options = "o:f:m:s:v";

static struct option long_options[] = {
    // general options
    {"output-dir",              required_argument,	     0,	         'o'},
    {"forget-param",            required_argument,       0,          'f'},
    {"frag-len-mean",		    required_argument,       0,          'm'},
    {"frag-len-std-dev",	    required_argument,       0,          's'},
    {"visualize",       	    no_argument,             0,          'v'},
    {"upper-quartile-norm",     no_argument,             0,          OPT_UPPER_QUARTILE_NORM},
    {"calc-covar",              no_argument,             0,          OPT_CALC_COVAR},
    {"no-bias-correct",         no_argument,             0,          OPT_NO_BIAS_CORRECT},
    {"in-between",              no_argument,             0,          OPT_IN_BETWEEN},
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

/**
 * Parse an float out of optarg and enforce that is between 'lower' and 'upper';
 * if it is not, output the given error message and
 * exit with an error and a usage message.
 */
float parseFloat(float lower, float upper, const char *errmsg, void (*print_usage)()) {
  float l;
  l = (float)atof(optarg);
  
  if (l < lower) {
    cerr << errmsg << endl;
    print_usage();
    exit(1);
  }
  
  if (l > upper)
    {
      cerr << errmsg << endl;
      print_usage();
      exit(1);
    }
  
  return l;
  
  cerr << errmsg << endl;
  print_usage();
  exit(1);
  return -1;
}

/**
 * Parse the options and set global variables appropriately.
 */
int parse_options(int argc, char** argv)
{
    int option_index = 0;
    int next_option;
	
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
        switch (next_option) {
			case -1:     /* Done with options. */
				break;
	                case 'f':
			  ff_param = parseFloat(0.5, 1, "-f/--forget-param arg must be between 0.5 and 1", print_usage);
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
            case 'v':
                vis = true;
                break;
			case OPT_UPPER_QUARTILE_NORM:
                upper_quart_norm = true;
                break;
            case OPT_CALC_COVAR:
                calc_covar = true;
                break;
            case OPT_NO_BIAS_CORRECT:
                bias_correct = false;
                break;
	case OPT_IN_BETWEEN:
	  in_between = true;
	  break;
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
    return 0;
}


/**
 * This function handles the probabilistic assignment of multi-mapped reads.  The marginal likelihoods are calculated for each mapping,
 * and the mass of the fragment is divided based on the normalized marginals to update the model parameters.
 * @param mass_n double specifying the logged fragment mass
 * @param frag_p pointer to the fragment to probabilistically assign
 * @param fld pointer to the fragment length distribution
 * @param bias_table pointer to the bias parameter table
 * @param mismatch_table pointer to the erorr model table
 * @param trans_table pointer to the transcript table
 */
void process_fragment(double mass_n, Fragment* frag_p, FLD* fld, BiasBoss* bias_table, MismatchTable* mismatch_table, TranscriptTable* trans_table)
{
    Fragment& frag = *frag_p;
    
    assert(frag.num_maps());
    
    frag.maps()[0]->mapped_trans->incr_bundle_counts();
    
    if (frag.num_maps()==1)
    // maps to a single location
    {
        const FragMap& m = *frag.maps()[0];
        Transcript* t  = m.mapped_trans;
        if (!t)
        {
            cerr << "ERROR: Transcript " << m.name << " not found in reference fasta.";
            exit(1);
        }
        
        //update parameters
        t->add_mass(0, mass_n);
        if (m.pair_status() == PAIRED)
                fld->add_val(m.length(), mass_n);
        
        if (bias_table)
            bias_table->update_observed(m, mass_n);
        
        mismatch_table->update(m, mass_n);
       
        delete frag_p;
        return;
    }
    
    // maps to multiple locations
    
    // calculate marginal likelihoods
    vector<double> likelihoods(frag.num_maps());
    double total_likelihood = 0.0;
    for(size_t i = 0; i < frag.num_maps(); ++i)
    {
        const FragMap& m = *frag.maps()[i];
        Transcript* t = m.mapped_trans;
        likelihoods[i] = t->log_likelihood(m);
        total_likelihood = (i) ? log_sum(total_likelihood, likelihoods[i]):likelihoods[i];
    }
    if (sexp(total_likelihood) == 0)
    {
        delete frag_p;
        return;
    }

    // normalize marginal likelihoods
    TransID bundle_rep = trans_table->get_trans_rep(frag.maps()[0]->trans_id);
    for(size_t i = 0; i < frag.num_maps(); ++i)
    {
        const FragMap& m = *frag.maps()[i];
        Transcript* t  = m.mapped_trans;
        double p = likelihoods[i]-total_likelihood;
        double mass_t = mass_n + p;
        
        if (i > 0)
        {
            TransID m_rep = trans_table->get_trans_rep(m.trans_id);
            if (m_rep != bundle_rep)
                bundle_rep = trans_table->merge_bundles(bundle_rep, m_rep);
        }
        
        if (sexp(p) == 0)
            continue;
        
        // update parameters
        t->add_mass(p, mass_n);
        
        if (m.pair_status() == PAIRED)
            fld->add_val(m.length(), mass_t);
 
        if (bias_table)
            bias_table->update_observed(m, mass_t);
        mismatch_table->update(m, mass_t);
        
        if (calc_covar)
        {
            for(size_t j = i+1; j < frag.num_maps(); ++j)
            {
                const FragMap& m2 = *frag.maps()[j];
                double p2 = likelihoods[j]-total_likelihood;
                if (sexp(p2) == 0)
                    continue;
                
                double covar = 2*mass_n + p + p2;
                trans_table -> update_covar(m.trans_id, m2.trans_id, covar); 
            }
        }
    }
        
    delete frag_p;
}

/**
 * This is the driver function for the main processing thread.  This function keeps track of the current fragment mass and sends fragments to be
 * processed once they are passed by the parsing thread.
 * @param map_parser the parsing object
 * @param trans_table pointer to the transcript table
 * @param fld pointer to the fragment length distribution
 * @param bias_table pointer to the bias parameter table
 * @param mismatch_table pointer to the erorr model table
 * @return the total number of fragments processed
 */
size_t threaded_calc_abundances(ThreadedMapParser& map_parser, TranscriptTable* trans_table, FLD* fld, BiasBoss* bias_table, MismatchTable* mismatch_table)
{
    ParseThreadSafety ts;
    ts.proc_lk.lock();
    ts.parse_lk.lock();
    boost::thread parse(&ThreadedMapParser::threaded_parse, &map_parser, &ts, trans_table);
    
    size_t n = 0;
    size_t fake_n = 0;
    double mass_n = 0;
    
    // For outputting rhos on log scale
    ofstream running_expr_file((output_dir + "/ce_running.expr").c_str());
    trans_table->output_header(running_expr_file);    
    int m = 0;
    int i = 1;

    // Geometric distribution used for 'in-between' tests
    srand ( time(NULL) );
    boost::math::geometric geom(.00002);
    cout.precision(20);
    while(true)
    {
        ts.proc_lk.lock();
        Fragment* frag = ts.next_frag;
        
        ts.parse_lk.unlock();
        if (!frag)
        {
            ts.proc_lk.unlock();
            break;
        }
        
        n++;

        //Output current values on log scale
        if (n == i * pow(10.,m))
        {
            running_expr_file << n << '\t';  
            trans_table->output_current(running_expr_file);
            running_expr_file.flush();
            m += (i==9);
            i = (i==9) ? 1 : (i+1);
        }
        
        // Output progress
        if (n == 1 || n % 10000 == 0)
        {
            if (vis)
            {
                cout << "0 " << fld->to_string() << '\n';
                cout << "1 " << mismatch_table->to_string() << '\n';
                cout << "2 " << bias_table->to_string() << '\n';
                cout << "3 " << n << '\n';
            }
            else
            {
                cout << n << '\t' << scientific << mass_n << '\t' << trans_table->covar_size() << '\t' << trans_table->num_bundles() <<'\n';
            }
        }
        
        
        // Update mass_n based on forgetting factor
        int k = (in_between) ? quantile(geom, rand()/double(RAND_MAX)) : 1;
        for (int j = 0; j < k; ++j)
        {
            if (++fake_n > 1)
            {
                mass_n += ff_param*log(fake_n-1) - log(pow(fake_n,ff_param) - 1);
            }
        }
        assert(!isnan(mass_n));
        
        process_fragment(mass_n, frag, fld, bias_table, mismatch_table, trans_table);
    }

    if (vis)
        cout << "99\n";
    else
        cout << "END: " << n << "\n";
        
    running_expr_file << n << '\t';
    trans_table->output_current(running_expr_file);
    running_expr_file.close();
    
    return n;
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
        
    FLD fld(fld_alpha, def_fl_max, def_fl_mean, def_fl_stddev);
    BiasBoss* bias_table = (bias_correct) ? new BiasBoss(bias_alpha):NULL;
    MismatchTable mismatch_table(mm_alpha);
    ThreadedMapParser map_parser(sam_hits_file_name);
    TranscriptTable trans_table(trans_fasta_file_name, expr_alpha, &fld, bias_table, &mismatch_table);

    if (bias_table)
        boost::thread bias_update(&TranscriptTable::threaded_bias_update, &trans_table);
    size_t tot_counts = threaded_calc_abundances(map_parser, &trans_table, &fld, bias_table, &mismatch_table);
    
    running = false;
    
    //mismatch_table.output(output_dir);
    ofstream fld_out((output_dir + "/fld.pdf").c_str());
    fld_out << fld.to_string() << '\n';
    fld_out.close();

    trans_table.output_bundles(output_dir);
    trans_table.output_expression(output_dir, tot_counts);
    return 0;
}

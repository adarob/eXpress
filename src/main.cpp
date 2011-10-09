//
//  main.cpp
//  express
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

//TODO: Use Bundle objects
//TODO: Use SAM/BAM indices
//TODO: Switch to compressed_matrix for covars
//TODO: Add additional EM step
//TODO: Speed up bias initialization (using non-logged values)

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "main.h"
#include "bundles.h"
#include "transcripts.h"
#include "fld.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include "mapparser.h"

#ifndef WIN32
    #include "update_check.h"
#endif

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

// the forgetting factor parameter controls the growth of the fragment mass
double ff_param = 0.9;

// the burn-in parameter determines how many reads are required before the 
// error and bias models are applied to probabilistic assignment 
size_t burn_in = 100000;

size_t stop_at = 0;

// file location parameters
string output_dir = ".";
string fasta_file_name = "";
string in_map_file_name = "";
string out_map_file_name = "";

// intial pseudo-count parameters (non-logged)
double expr_alpha = .001;
double fld_alpha = 1;
double bias_alpha = 1;
double mm_alpha = 1;

// fragment length parameters
int def_fl_max = 800;
int def_fl_mean = 200;
int def_fl_stddev = 60;

// option parameters
bool error_model = true;
bool bias_correct = true;
bool calc_covar = false;
bool vis = false;
bool output_alignments = false;
bool output_running = false;

// directional parameters
enum Direction{ BOTH, FR, RF };
Direction direction = BOTH;

bool running = true;

// used for multiple rounds of EM
enum Iteration{ FIRST, LAST, ONLY };
Iteration iteration = ONLY;

bool parse_options(int ac, char ** av)
{
    po::options_description generic("Allowed options");
    generic.add_options()
    ("help,h", "produce help message")
    ("output-dir,o", po::value<string>(&output_dir)->default_value("."), "write all output files to this directory")
    ("frag-len-mean,m", po::value<int>(&def_fl_mean)->default_value(200), "prior estimate for average fragment length")
    ("frag-len-stddev,s", po::value<int>(&def_fl_stddev)->default_value(60), "prior estimate for fragment length std deviation")
    ("output-alignments", "output alignments (sam/bam) with probabilistic assignments")
    ("fr-stranded", "accept only forward->reverse alignments (second-stranded protocols)")
    ("rf-stranded", "accept only reverse->forward alignments (first-stranded protocols)")
    ("calc-covar", "calculate and output covariance matrix")
    ("no-update-check", "disables automatic check for update via web")
    ;
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
    ("no-bias-correct","")
    ("no-error-model","")
    ("single-round", "")
    ("visualizer-output,v","")
    ("output-running", "")
    ("forget-param,f", po::value<double>(&ff_param)->default_value(0.9),"")
    ("stop-at", po::value<size_t>(&stop_at)->default_value(0),"")
    ("sam-file", po::value<string>(&in_map_file_name)->default_value(""),"")
    ("fasta-file", po::value<string>(&fasta_file_name)->default_value(""),"")
    ;
    
    po::positional_options_description positional;
    positional.add("fasta-file",1).add("sam-file",1);
    
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    
    bool error = false;
    po::variables_map vm;
    try 
    {
        po::store(po::command_line_parser(ac, av).options(cmdline_options).positional(positional).run(), vm);
    }
    catch (po::error& e)
    {
        cerr << "Command-Line Argument Error: "<< e.what() << endl;
        error = true;
    }
    po::notify(vm);
    
    if (fasta_file_name == "")
    {
        cerr << "Command-Line Argument Error: target sequence fasta file required\n\n";
        error = true;
    }
    
    if (error || vm.count("help")) 
    {
        cerr << "express v" << PACKAGE_VERSION << endl;
        cerr << "-----------------------------\n"; 
        cerr << "File Usage:  express [options] <target_seqs.fa> <hits.(sam/bam)>\n";
        cerr << "Piped Usage: bowtie [options] -S <index> <reads.fq> | express [options] <target_seqs.fa>\n";
        cerr << "Required arguments:\n";
        cerr << " <target_seqs.fa>                     target sequence file in fasta format\n";
        cerr << " <hits.(sam/bam)>                     read alignment file in SAM or BAM format\n";
        
        cerr << generic;
        
        return 1;
    }
    
    if (vm.count("fr-stranded"))
    {
        direction = FR;
    }
    
    if (vm.count("rf-stranded"))
    {
        if (direction != BOTH)
        {
            cerr << "fr-stranded and rf-stranded flags cannot both be specified in the same run.\n";
            return 1;
        }
        direction = RF;
    }
    
    calc_covar = vm.count("calc-covar");
    bias_correct = !(vm.count("no-bias-correct"));
    error_model = !(vm.count("no-error-model"));
    vis = vm.count("visualizer-output");
    output_running = vm.count("output-running");
    
    if (!vm.count("single-round") && in_map_file_name != "")
        iteration = FIRST;

    if (vm.count("output-alignments"))
    {
        out_map_file_name = output_dir + "hits.prob";
    }
    
#ifndef WIN32
    if (!vm.count("no-update-check"))
        check_version(PACKAGE_VERSION);
#endif
    
    return 0;
}


/**
 * This function handles the probabilistic assignment of multi-mapped reads.  The marginal likelihoods are calculated for each mapping,
 * and the mass of the fragment is divided based on the normalized marginals to update the model parameters.
 * @param mass_n double specifying the logged fragment mass
 * @param frag_p pointer to the fragment to probabilistically assign
 * @param trans_table pointer to the transcript table
 * @param globs a pointer to the struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
 */
void process_fragment(double mass_n, Fragment* frag_p, TranscriptTable* trans_table, const Globals& globs)
{
    Fragment& frag = *frag_p;
    
    assert(frag.num_hits());
    
    Bundle* bundle = frag.hits()[0]->mapped_trans->bundle();
    if (iteration != FIRST)
        bundle->incr_counts();
    
    if (frag.num_hits()==1)
    // hits to a single location
    {
        FragHit& m = *frag.hits()[0];
        Transcript* t  = m.mapped_trans;
        
        m.probability = 1.0;
        
        //update parameters
        if (iteration == FIRST || iteration == ONLY)
        {
            t->add_mass(0, mass_n);
            t->incr_uniq_counts();
            if (m.pair_status() == PAIRED)
                    (globs.fld)->add_val(m.length(), mass_n);
            
            if (globs.bias_table)
                (globs.bias_table)->update_observed(m, mass_n - t->mass() + log((double)t->unbiased_effective_length()));
            if (globs.mismatch_table)
                (globs.mismatch_table)->update(m, mass_n);
        }
        else
        {
            t->add_prob_count(1.0);
        }
        return;
    }
    
    // hits to multiple locations
    
    // calculate marginal likelihoods
    vector<double> likelihoods(frag.num_hits());
    double total_likelihood = 0.0;
    for(size_t i = 0; i < frag.num_hits(); ++i)
    {
        const FragHit& m = *frag.hits()[i];
        Transcript* t = m.mapped_trans;
        likelihoods[i] = t->log_likelihood(m);
        total_likelihood = (i) ? log_sum(total_likelihood, likelihoods[i]):likelihoods[i];
    }
    if (sexp(total_likelihood) == 0)
    {
        return;
    }
    
    // merge bundles
    for(size_t i = 1; i < frag.num_hits(); ++i)
    {
        bundle = trans_table->merge_bundles(bundle, frag.hits()[i]->mapped_trans->bundle());
    }

    // normalize marginal likelihoods
    for(size_t i = 0; i < frag.num_hits(); ++i)
    {
        FragHit& m = *frag.hits()[i];
        Transcript* t  = m.mapped_trans;
        double p = likelihoods[i]-total_likelihood;
        double mass_t = mass_n + p;
               
        m.probability = sexp(p);
        if (m.probability == 0)
            continue;
        assert(!isinf(m.probability));
        
        // update parameters
        if (iteration == FIRST || iteration == ONLY)
        {

            t->add_mass(p, mass_n);
            
            if (m.pair_status() == PAIRED)
                (globs.fld)->add_val(m.length(), mass_t);
     
            if (globs.bias_table)
                (globs.bias_table)->update_observed(m, mass_t - t->mass() + log((double)t->unbiased_effective_length()));
            if (globs.mismatch_table)
                (globs.mismatch_table)->update(m, mass_t);
        }
        else
        {
            t->add_prob_count(m.probability);
        }
        
        if (calc_covar && (iteration == ONLY || iteration == LAST))
        {
            for(size_t j = i+1; j < frag.num_hits(); ++j)
            {
                const FragHit& m2 = *frag.hits()[j];
                double p2 = likelihoods[j]-total_likelihood;
                if (sexp(p2) == 0)
                    continue;
                
                double covar = p + p2;
                if (iteration == ONLY)
                    covar += 2*mass_n;
                
                trans_table->update_covar(m.trans_id, m2.trans_id, covar); 
            }
        }
    }
}

/**
 * This is the driver function for the main processing thread.  This function keeps track of the current fragment mass and sends fragments to be
 * processed once they are passed by the parsing thread.
 * @param map_parser the parsing object
 * @param trans_table pointer to the transcript table
 * @param globs a struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
 * @return the total number of fragments processed
 */
size_t threaded_calc_abundances(ThreadedMapParser& map_parser, TranscriptTable* trans_table, Globals& globs)
{
    cout << "Processing input fragment alignments...\n";
    running = true;
    
    ParseThreadSafety ts;
    boost::thread parse(&ThreadedMapParser::threaded_parse, &map_parser, &ts, trans_table);
    boost::thread* bias_update = NULL;
    
    size_t n = 0;
    double mass_n = 0;
    Fragment* frag;
    cout << setiosflags(ios::left);
    
    while(!stop_at || n < stop_at)
    {
        {
            boost::unique_lock<boost::mutex> lock(ts.mut);
            while(!ts.frag_clean)
            {
                ts.cond.wait(lock);
            }
            
            frag = ts.next_frag;
            
            ts.frag_clean = false;
            ts.cond.notify_one();
        }
        
        if (!frag)
        {
            break;
        }
        
        n++;
        
        if (n == burn_in)
        {
            if (globs.bias_table)
                bias_update = new boost::thread(&TranscriptTable::threaded_bias_update, trans_table);
            if (globs.mismatch_table)
                (globs.mismatch_table)->activate();
        }
        
        // Output progress
        if (n == 1 || n % 100000 == 0)
        {
            if (output_running && (iteration == ONLY || iteration == LAST))
            {
                char buff[500];
                sprintf(buff, "%s/x_" SIZE_T_FMT "", output_dir.c_str(), n);
                string dir(buff);
                try { fs::create_directories(dir); }
                catch (fs::filesystem_error& e)
                {
                    cerr << e.what() << endl;
                    exit(1);
                }
                trans_table->output_results(dir, n, false, iteration==LAST);
                ofstream paramfile((dir + "/params.xprs").c_str());
                (globs.fld)->append_output(paramfile);
                if (globs.mismatch_table)
                    (globs.mismatch_table)->append_output(paramfile);
                if (globs.bias_table)
                    (globs.bias_table)->append_output(paramfile);
                paramfile.close();
            }
            if (vis)
            {
                cout << "0 " << (globs.fld)->to_string() << '\n';
                cout << "1 " << (globs.mismatch_table)->to_string() << '\n';
                cout << "2 " << (globs.bias_table)->to_string() << '\n';
                cout << "3 " << n << '\n';
            }
            else
            {
                cout << "Fragments Processed: " << setw(9) << n << "\t Number of Bundles: "<< trans_table->num_bundles() << endl;
            }
        }
        
        
        if (n > 1)
        {
            mass_n += ff_param*log((double)n-1) - log(pow(n,ff_param) - 1);
        }
        
        process_fragment(mass_n, frag, trans_table, globs);
    }
    
    {
        boost::unique_lock<boost::mutex> lock(ts.mut);
        running = false;
        ts.frag_clean = false;
        ts.cond.notify_one();
    }
    
    if (vis)
        cout << "99\n";
    else
        cout << "COMPLETED: Processed " << n << " fragments, targets are in " << trans_table->num_bundles() << " bundles\n";
    
    parse.join();
    if (bias_update)
    {
        bias_update->join();
        delete bias_update;
    }
    return n;
}


int main (int argc, char ** argv)
{      
    int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
    
    if (output_dir != ".")
    {
        try { fs::create_directories(output_dir); }
        catch (fs::filesystem_error& e)
        {
            cerr << e.what() << endl;; 
        }
    }
    if (!fs::exists(output_dir))
    {
        cerr << "ERROR: cannot create directory " << output_dir << ".\n";
        exit(1);
    }
    
    Globals globs;
    globs.fld = new FLD(fld_alpha, def_fl_max, def_fl_mean, def_fl_stddev);
    globs.bias_table = (bias_correct) ? new BiasBoss(bias_alpha):NULL;
    globs.mismatch_table = (error_model) ? new MismatchTable(mm_alpha):NULL;
    ThreadedMapParser* map_parser = new ThreadedMapParser(in_map_file_name, (iteration==ONLY) ? out_map_file_name:"");
    TranscriptTable trans_table(fasta_file_name, map_parser->trans_index(), expr_alpha, &globs);

    double num_trans = (double)map_parser->trans_index().size();
    if (calc_covar && (double)SIZE_MAX < num_trans*(num_trans+1))
    {
        cerr << "Warning: Your system is unable to represent large enough values for efficiently hashing transcript pairs.  Covariance calculation will be disabled.\n";
        calc_covar = false;
    }
    
    size_t tot_counts = threaded_calc_abundances(*map_parser, &trans_table, globs);
    delete map_parser;
    
    if (iteration == FIRST)
    {
        cout << "\nRe-estimating counts with second-round of EM...\n";
        iteration = LAST;
        map_parser = new ThreadedMapParser(in_map_file_name, out_map_file_name);
        tot_counts = threaded_calc_abundances(*map_parser, &trans_table, globs);
    }
    
	cout << "Writing results to file...\n";
    
    trans_table.output_results(output_dir, tot_counts, calc_covar, iteration!=ONLY);
    ofstream paramfile((output_dir + "/params.xprs").c_str());
    (globs.fld)->append_output(paramfile);
    if (globs.mismatch_table)
        (globs.mismatch_table)->append_output(paramfile);
    if (globs.bias_table)
        (globs.bias_table)->append_output(paramfile);
    paramfile.close();
    
    cout << "Done\n";
    
    return 0;
}

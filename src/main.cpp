//
//  main.cpp
//  express
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "main.h"
#include "bundles.h"
#include "targets.h"
#include "fld.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include "mapparser.h"
#include "threadsafety.h"
#include "robertsfilter.h"

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
//size_t burn_in = 100;
size_t burn_out = 5000000;
bool burned_out = false;

size_t stop_at = 0;

// file location parameters
string output_dir = ".";
string fasta_file_name = "";
string in_map_file_name = "";
string out_map_file_name = "";

// intial pseudo-count parameters (non-logged)
double expr_alpha = .1;
double fld_alpha = 1;
double bias_alpha = 1;
double mm_alpha = 1;

// fragment length parameters
int def_fl_max = 800;
int def_fl_mean = 200;
int def_fl_stddev = 80;

// option parameters
bool edit_detect = false;
bool error_model = true;
bool bias_correct = true;
bool calc_covar = false;
bool output_alignments = false;
bool output_running_rounds = false;
bool output_running_reads = false;
size_t num_threads = 2;

// directional parameters
enum Direction{ BOTH, FR, RF };
Direction direction = BOTH;

bool running = true;

// used for multiple rounds of EM
bool first_round = true;
bool last_round = true;
bool batch_mode = false;
bool online_additional = false;
bool both = false;
size_t remaining_rounds = 0;

typedef boost::unordered_map<string, double> AlphaMap;
AlphaMap* expr_alpha_map = NULL;

AlphaMap* parse_priors(string in_file)
{
    ifstream ifs(in_file.c_str());
    if(!ifs.is_open())
    {
        cerr << "ERROR: Unable to open input priors file '" << in_file << "'.\n" ; 
        exit(1);
    }
    AlphaMap* alphas = new AlphaMap();
    
    string line;
    
    while(ifs.good())
    {
        getline(ifs,line);
        
        size_t idx = line.find_first_of("\t ");
        if (idx!=string::npos)
        {
            string name = line.substr(0,idx);
            string val = line.substr(idx+1);
            (*alphas)[name] = atof(val.c_str());
        }
    }
    return alphas;
};

bool parse_options(int ac, char ** av)
{
    po::options_description generic("Allowed options");
    generic.add_options()
    ("help,h", "produce help message")
    ("output-dir,o", po::value<string>(&output_dir)->default_value("."), "write all output files to this directory")
    ("num-threads,p", po::value<size_t>(&num_threads)->default_value(2), "number of threads (>= 2)")
    ("frag-len-mean,m", po::value<int>(&def_fl_mean)->default_value(def_fl_mean), "prior estimate for average fragment length")
    ("frag-len-stddev,s", po::value<int>(&def_fl_stddev)->default_value(def_fl_stddev), "prior estimate for fragment length std deviation")
    ("additional-batch,B", po::value<size_t>(&remaining_rounds)->default_value(remaining_rounds), "number of additional batch EM rounds after initial online round")
    ("additional-online,O", po::value<size_t>(&remaining_rounds)->default_value(remaining_rounds), "number of additional online EM rounds after initial online round")
    ("output-align-prob", "output alignments (sam/bam) with probabilistic assignments")
    ("output-align-samp", "output alignments (sam/bam) with sampled assignments")
    ("fr-stranded", "accept only forward->reverse alignments (second-stranded protocols)")
    ("rf-stranded", "accept only reverse->forward alignments (first-stranded protocols)")
    ("calc-covar", "calculate and output covariance matrix")
    ("no-update-check", "disables automatic check for update via web")
    ;
    
    string prior_file = "";
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
    ("edit-detect","")
    ("no-bias-correct","")
    ("no-error-model","")
    ("single-round", "")
    ("output-running-rounds", "")    
    ("output-running-reads", "")
    ("batch-mode","")
    ("both","")
    ("burn-out", po::value<size_t>(&burn_out)->default_value(burn_out), "")
    ("prior-params", po::value<string>(&prior_file)->default_value(""), "")
    ("forget-param,f", po::value<double>(&ff_param)->default_value(ff_param),"")
    ("expr-alpha", po::value<double>(&expr_alpha)->default_value(expr_alpha),"")
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
    
    if (ff_param > 1.0 || ff_param < 0.5)
    {
        cerr << "Command-Line Argument Error: forget-param/f option must be between 0.5 and 1.0\n\n";
        error= true;
    }
        
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

    edit_detect = vm.count("edit-detect");
    calc_covar = vm.count("calc-covar");
    bias_correct = !(vm.count("no-bias-correct"));
    error_model = !(vm.count("no-error-model"));
    output_running_rounds = vm.count("output-running-rounds");
    output_running_reads = vm.count("output-running-reads");
    batch_mode = vm.count("batch-mode");
    online_additional = vm.count("additional-online");
    both = vm.count("both");
    
    if (num_threads < 2)
        num_threads = 0;
    num_threads -= 2;
    if (num_threads > 0)
        num_threads -= edit_detect;
    
    if (remaining_rounds > 0 && in_map_file_name != "")
        last_round = false;

    if (vm.count("output-align-prob"))
    {
        out_map_file_name = output_dir + "/hits.prob";
    }
    
    if (vm.count("output-align-samp"))
    {
        out_map_file_name = output_dir + "/hits.samp";
    }
    
    if (prior_file != "")
    {
        expr_alpha_map = parse_priors(prior_file);
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
 * @param frag_p pointer to the fragment to probabilistically assign
 * @param globs a pointer to the struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
 */
void process_fragment(Fragment* frag_p, Globals& globs)
{
    Fragment& frag = *frag_p;
    frag.sort_hits();
    double mass_n = frag.mass();
    
    assert(frag.num_hits());
        
    // calculate marginal likelihoods
    vector<double> likelihoods(frag.num_hits(),0);
    vector<double> masses(frag.num_hits(),0);
    vector<double> variances(frag.num_hits(), 0);
    double total_likelihood = HUGE_VAL;
    double total_mass = HUGE_VAL;
    double total_variance = HUGE_VAL;
    size_t num_solvable = 0;
    
    
    boost::unordered_set<Target*> targ_set;

    if (frag.num_hits()>1)
    {
        for(size_t i = 0; i < frag.num_hits(); ++i)
        {
            
            const FragHit& m = *frag.hits()[i];
            Target* t = m.mapped_targ;
            if (targ_set.count(t) == 0)
            {
                t->lock();
                targ_set.insert(t);
            }
            likelihoods[i] = t->log_likelihood(m, first_round);
            masses[i] = t->mass();
            variances[i] = t->mass_var();
            total_likelihood = log_sum(total_likelihood, likelihoods[i]);
            total_mass = log_sum(total_mass, masses[i]);
            total_variance = log_sum(total_variance, variances[i]);
            num_solvable += t->solvable();
        }
    }
    else
    {
        Target* t = frag.hits()[0]->mapped_targ;
        t->lock();
        targ_set.insert(t);
        total_likelihood = likelihoods[0];
    }
    
    assert(!islzero(total_likelihood));
    
    // merge bundles
    Bundle* bundle = frag.hits()[0]->mapped_targ->bundle();
    if (first_round)
    {
        bundle->incr_counts();
    }
    
    // normalize marginal likelihoods
    for(size_t i = 0; i < frag.num_hits(); ++i)
    {
        FragHit& m = *frag.hits()[i];
        Target* t  = m.mapped_targ;
 
        bundle = globs.targ_table->merge_bundles(bundle, t->bundle());
        
        double p = likelihoods[i]-total_likelihood;
        double mass_t = mass_n + p;
        
        double v = HUGE_VAL;
        if (frag.num_hits() > 1)
            v = log_sum(variances[i] - 2*total_mass, total_variance + 2*masses[i] - 4*total_mass);
        
        assert(!isnan(v));
        assert(!(isnan(mass_t)||isinf(mass_t)));
        
        m.probability = sexp(p);
        
        assert(!isinf(m.probability));
        
        t->add_mass(p, v, mass_n);
        
        // update parameters
        if (first_round)
        {
            t->incr_counts(frag.num_hits()==1);
            if (!t->solvable() && num_solvable == frag.num_hits()-1)
            {
                t->solvable(true);
            }
            if ((!burned_out || edit_detect) && globs.mismatch_table)
                    (globs.mismatch_table)->update(m, p, mass_n);
            if (!burned_out)
            {
                if (m.pair_status() == PAIRED)
                    (globs.fld)->add_val(m.length(), mass_t);
                if (globs.bias_table)
                    (globs.bias_table)->update_observed(m, mass_t);
            }
        }
        
        if (calc_covar && (last_round || online_additional))
        {
            for(size_t j = i+1; j < frag.num_hits(); ++j)
            {
                const FragHit& m2 = *frag.hits()[j];
                double p2 = likelihoods[j]-total_likelihood;
                if (sexp(p2) == 0)
                    continue;
                
                double covar = p + p2;
                if ((first_round && last_round) || online_additional)
                    covar += 2*mass_n;
                
                globs.targ_table->update_covar(m.targ_id, m2.targ_id, covar); 
            }
        }
    }
    
    foreach (Target* t, targ_set)
    {
        t->unlock();
    }
}

void proc_thread(ParseThreadSafety* pts, Globals* globs)
{
    while (true)
    {
        Fragment* frag = pts->proc_on.pop();
        if (!frag)
            break;
        process_fragment(frag, *globs);
        pts->proc_out.push(frag);
    }
}

/**
 * This is the driver function for the main processing thread.  This function keeps track of the current fragment mass and sends fragments to be
 * processed once they are passed by the parsing thread.
 * @param map_parser the parsing object
 * @param targ_table pointer to the target table
 * @param globs a struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
 * @return the total number of fragments processed
 */
size_t threaded_calc_abundances(ThreadedMapParser& map_parser, TargetTable* targ_table, Globals& globs)
{
    cout << "Processing input fragment alignments...\n";
    
    vector<boost::thread*> thread_pool;
    boost::thread* bias_update = NULL;
    
    RobertsFilter frags_seen;
    
    size_t n = 1;
    size_t num_frags = 0;
    double mass_n = 0;
    cout << setiosflags(ios::left);
    
    // For log-scale output
    size_t i = 1;
    size_t j = 6;
    
    
    Fragment* frag;
    
    while (true)
    {
        boost::mutex bu_mut;
        running = true;
        ParseThreadSafety pts(max((int)num_threads,10));
        boost::thread parse(&ThreadedMapParser::threaded_parse, &map_parser, &pts, targ_table, stop_at);
        
        if (burned_out)
        {
            for(size_t i=0; i < thread_pool.size(); i++)
            {
                thread_pool[i] = new boost::thread(proc_thread, &pts, &globs);
            }
        }
        
        while(true)
        {
            frag = pts.proc_in.pop();
            if (frag && frags_seen.test_and_push(frag->name()))
            {
                cerr << "ERROR: Alignments are not properly sorted. Read '" << frag->name() << "' has alignments which are non-consecutive.\n" ; 
                exit(1);  
            }
            if (num_threads && burned_out)
            {
                if (!frag)
                {
                    for(size_t i = 0; i < thread_pool.size(); ++i)
                    {
                        pts.proc_on.push(NULL);
                    }
                    break;
                }
                frag->mass(mass_n);
                pts.proc_on.push(frag);
            }
            else
            {
                if (!frag)
                {
                    break;
                }
                {
                    frag->mass(mass_n);
                    boost::unique_lock<boost::mutex> lock(bu_mut);
                    process_fragment(frag, globs);
                    pts.proc_out.push(frag);
                }
            }
                    
            if (n == burn_in)
            {
                bias_update = new boost::thread(&TargetTable::threaded_bias_update, targ_table, &bu_mut);
                if (globs.mismatch_table)
                    (globs.mismatch_table)->activate();
            }
            if (n == burn_out)
            {
                (globs.mismatch_table)->fix();
                burned_out = true;
                thread_pool = vector<boost::thread*>(num_threads);
                for(size_t i=0; i < thread_pool.size(); i++)
                {
                    thread_pool[i] = new boost::thread(proc_thread, &pts, &globs);
                }
            }
            
            if (output_running_reads && n == i*pow(10.,(double)j))
            {
                boost::unique_lock<boost::mutex> lock(bu_mut);
                char buff[500];
                sprintf(buff, "%s/x_" SIZE_T_FMT "", output_dir.c_str(), n);
                cout << "Writing results to " << buff << endl;
                string dir(buff);
                try { fs::create_directories(dir); }
                catch (fs::filesystem_error& e)
                {
                    cerr << e.what() << endl;
                    exit(1);
                }
                targ_table->output_results(dir, n);
                
                ofstream paramfile((dir + "/params.xprs").c_str());
                (globs.fld)->append_output(paramfile);
                if (globs.mismatch_table)
                    (globs.mismatch_table)->append_output(paramfile);
                if (globs.bias_table)
                    (globs.bias_table)->append_output(paramfile);
                paramfile.close();

                if (i++ == 9)
                {
                    i = 1;
                    j++;
                }
            }
            
            num_frags++;
            
            // Output progress
            if (num_frags % 1000000 == 0)
            {
                cout << "Fragments Processed: " << setw(9) << num_frags << "\t Number of Bundles: "<< targ_table->num_bundles() << endl;
            }
            n++;
            mass_n += ff_param*log((double)n-1) - log(pow(n,ff_param) - 1);
        }
    
        running = false;

        parse.join();
        foreach(boost::thread* t, thread_pool)
        {
            t->join();
        }
        
        if (bias_update)
        {
            bias_update->join();
            delete bias_update;
            bias_update = NULL;
        }
        
        if (online_additional && remaining_rounds--)
        {
            if (output_running_rounds)
            {
                char buff[500];
                sprintf(buff, "%s/x_" SIZE_T_FMT "", output_dir.c_str(), remaining_rounds);
                string dir(buff);
                try { fs::create_directories(dir); }
                catch (fs::filesystem_error& e)
                {
                    cerr << e.what() << endl;
                    exit(1);
                }
                targ_table->output_results(dir, num_frags);
            }
            cout << remaining_rounds << " remaining rounds." << endl;
            first_round = false;
            last_round = (remaining_rounds==0 && !both);
            map_parser.write_active(last_round);
            map_parser.reset_reader();
            num_frags = 0;
        }
        else
        {
            break;
        }
    }
        
 
    cout << "COMPLETED: Processed " << num_frags << " mapped fragments, targets are in " << targ_table->num_bundles() << " bundles\n";
    
    return num_frags;
}


int main (int argc, char ** argv)
{     
    srand ( time(NULL) );
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
    globs.mismatch_table = (error_model) ? new MismatchTable(mm_alpha):NULL;
    
    ThreadedMapParser map_parser(in_map_file_name, out_map_file_name, last_round);
    globs.bias_table = (bias_correct) ? new BiasBoss(bias_alpha):NULL;
    TargetTable targ_table(fasta_file_name, map_parser.targ_index(), map_parser.targ_lengths(), edit_detect, expr_alpha, expr_alpha_map, &globs);
    globs.targ_table = &targ_table;
    
    double num_targ = (double)map_parser.targ_index().size();
    
    if (calc_covar && (double)SSIZE_MAX < num_targ*(num_targ+1))
    {
        cerr << "Warning: Your system is unable to represent large enough values for efficiently hashing target pairs.  Covariance calculation will be disabled.\n";
        calc_covar = false;
    }
    
    size_t tot_counts = threaded_calc_abundances(map_parser, &targ_table, globs);
    
    if (both)
    {
        remaining_rounds = 1;
        online_additional = false;
    }
    
    targ_table.round_reset();
    ff_param = 1.0;
    first_round = false;
    while (!last_round)
    {
        if (output_running_rounds)
        {
            char buff[500];
            sprintf(buff, "%s/x_" SIZE_T_FMT "", output_dir.c_str(), remaining_rounds);
            string dir(buff);
            try { fs::create_directories(dir); }
            catch (fs::filesystem_error& e)
            {
                cerr << e.what() << endl;
                exit(1);
            }
            targ_table.output_results(dir, tot_counts);
        }
        remaining_rounds--;
        cout << "\nRe-estimating counts with additional round of EM (" << remaining_rounds << " remaining)...\n";
        last_round = (remaining_rounds == 0);
        map_parser.write_active(last_round);
        map_parser.reset_reader();
        tot_counts = threaded_calc_abundances(map_parser, &targ_table, globs);
        targ_table.round_reset();
    }
    
	cout << "Writing results to file...\n";
    
    targ_table.output_results(output_dir, tot_counts, calc_covar, edit_detect);
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

//
//  main.cpp
//  express
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

//TODO: Ouptut covariances
//TODO: Output confidence intervals

//#include <boost/math/distributions/geometric.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
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
#include "version.h"
#include "update_check.h"

using namespace std;
namespace po = boost::program_options;

// the forgetting factor parameter controls the growth of the fragment mass
double ff_param = 0.9;

string output_dir = ".";
string fasta_file_name = "";
string sam_file_name = "";


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
bool output_running = false;

//bool in_between = false;
enum Direction{ BOTH, FR, RF };
Direction direction = BOTH;

bool running = true;

bool parse_options(int ac, char ** av)
{
    po::options_description generic("Allowed options");
    generic.add_options()
    ("help,h", "produce help message")
    ("output-dir,o", po::value<string>(&output_dir)->default_value("."), "write all output files to this directory")
    ("frag-len-mean,m", po::value<int>(&def_fl_mean)->default_value(200), "prior estimate for average fragment length")
    ("frag-len-stddev,s", po::value<int>(&def_fl_stddev)->default_value(60), "prior estimate for fragment length std deviation")
    ("fr-stranded", "accept only forward->reverse alignments (second-stranded protocols)")
    ("rf-stranded", "accept only reverse->forward alignments (first-stranded protocols)")
    ("calc-covar", "calculate and output covariance matrix")
    ("no-update-check", "disables auotmatic check for update via web")
    ;
    
    po::options_description hidden("Hidden options");
    hidden.add_options()
    ("no-bias-correct","")
    ("no-error-model","")
      //    ("in-between","")
    ("visualizer-output,v","")
    ("output-running","")
    ("forget-param,f", po::value<double>(&ff_param)->default_value(0.9),"")
    ("sam-file", po::value<string>(&sam_file_name)->default_value(""),"")
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
    //  in_between = vm.count("in-between");
    vis = vm.count("visualizer-output");
    output_running = vm.count("output-running");
    
    if (!vm.count("no-update-check"))
    {
        check_version(PACKAGE_VERSION);
    }
    
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
        if (mismatch_table)
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
        if (mismatch_table)
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
    cout << "Processing input fragment alignments...\n";
    
    ParseThreadSafety ts;
    ts.proc_lk.lock();
    ts.parse_lk.lock();
    boost::thread parse(&ThreadedMapParser::threaded_parse, &map_parser, &ts, trans_table);
    
    size_t n = 0;
    size_t fake_n = 0;
    double mass_n = 0;
    
    // For outputting rhos on log scale
    int m = 0;
    int i = 1;
    ofstream running_expr_file;
    if (output_running)
    {
        running_expr_file.open((output_dir + "/running.xprs").c_str());
        trans_table->output_header(running_expr_file);    
    }
    
    // Geometric distribution used for 'in-between' tests
    //srand ( time(NULL) );
    //boost::math::geometric geom(.00002);
    cout << setiosflags(ios::left);
    
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
        if (output_running && n == i * pow(10.,m))
        {
            running_expr_file << n << '\t';  
            trans_table->output_current(running_expr_file);
            running_expr_file.flush();
            m += (i==9);
            i = (i==9) ? 1 : (i+1);
        }
        
        // Output progress
        if (n == 1 || n % 100000 == 0)
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
                cout << "Fragments Processed: " << setw(9) << n << "\t Number of Bundles: "<< trans_table->num_bundles() << endl;
            }
        }
        
        
        // Update mass_n based on forgetting factor
        //int k = (in_between) ? quantile(geom, rand()/double(RAND_MAX)) : 1;
        //for (int j = 0; j < k; ++j)
        //{
            if (++fake_n > 1)
            {
                mass_n += ff_param*log((double)fake_n-1) - log(pow(fake_n,ff_param) - 1);
            }
        //}
        assert(!isnan(mass_n));
        
        process_fragment(mass_n, frag, fld, bias_table, mismatch_table, trans_table);
    }

    if (vis)
        cout << "99\n";
    else
        cout << "COMPLETED: Processed " << n << " fragments, targets are in " << trans_table->num_bundles() << " bundles\n";
    
    if (output_running)
    {
        running_expr_file << n << '\t';
        trans_table->output_current(running_expr_file);
        running_expr_file.close();
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
        boost::filesystem::create_directories(output_dir);
    }
    if (!boost::filesystem::exists(output_dir))
    {
        cerr << "Error: cannot create directory " << output_dir << ".\n";
        exit(1);
    }
            
    FLD fld(fld_alpha, def_fl_max, def_fl_mean, def_fl_stddev);
    BiasBoss* bias_table = (bias_correct) ? new BiasBoss(bias_alpha):NULL;
    MismatchTable* mismatch_table = (error_model) ? new MismatchTable(mm_alpha):NULL;
    ThreadedMapParser map_parser(sam_file_name);
    TranscriptTable trans_table(fasta_file_name, expr_alpha, &fld, bias_table, mismatch_table);

    if (bias_table)
        boost::thread bias_update(&TranscriptTable::threaded_bias_update, &trans_table);

    size_t tot_counts = threaded_calc_abundances(map_parser, &trans_table, &fld, bias_table, mismatch_table);
	running = false;

	cout << "Writing results to file...\n";
    
    //mismatch_table->output(output_dir);
    //ofstream fld_out((output_dir + "/fld.out").c_str());
    //fld_out << fld.to_string() << '\n';
    //fld_out.close();
	//trans_table.output_bundles(output_dir);
    
	trans_table.output_expression(output_dir, tot_counts);
    return 0;
}

//
//  targets.cpp
//  express
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "main.h"
#include "targets.h"
#include "fld.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>

using namespace std;

Target::Target(const TargID id, const std::string& name, const std::string& seq, bool prob_seq, double alpha, const Globals* globs)
:   _globs(globs),
    _id(id),
    _name(name),
    _seq_f(seq, 0, prob_seq),
    _seq_r(_seq_f),
    _alpha(log(alpha)),
    _ret_params(&_curr_params),
    _uniq_counts(0),
    _tot_counts(0),
    _avg_bias(0),
    _solvable(false)
{ 
    if (globs->bias_table)
    {
        _start_bias = new std::vector<float>(seq.length(),0);
        _end_bias = new std::vector<float>(seq.length(),0);
    }
    else
    {
        _start_bias = NULL;
        _end_bias = NULL;
    }
    _cached_eff_len = est_effective_length(NULL, false);
    _curr_params.mass_var = _cached_eff_len + _alpha;
}

void Target::add_mass(double p, double v, double mass) 
{ 
    _curr_params.mass = log_sum(_curr_params.mass, p+mass);
    _curr_params.mass_var = log_sum(_curr_params.mass_var, p+mass*2);
    if (p != 0.0 || v != HUGE_VAL)
    {
        _curr_params.mass_var = log_sum(_curr_params.mass_var, v+2*mass);
        _curr_params.var_sum = log_sum(_curr_params.var_sum, v+mass);
        if (p!=HUGE_VAL)
            _curr_params.tot_ambig_mass = log_sum(_curr_params.tot_ambig_mass, mass);
    }
    
    (_globs->targ_table)->update_total_fpb(mass - _cached_eff_len);
    
}  

void Target::round_reset()
{
    _last_params = _curr_params;
    _curr_params = RoundParams();
    _ret_params = &_last_params;
}

double Target::rho() const
{
    double eff_len = cached_effective_length(false);
    if (eff_len == HUGE_VAL)
        return HUGE_VAL;
    
    return mass(true) - eff_len - (_globs->targ_table)->total_fpb();
}

double Target::mass(bool with_pseudo) const
{
    if (!with_pseudo)
        return _ret_params->mass;
    return log_sum(_ret_params->mass, _alpha+_cached_eff_len);
}

double Target::mass_var(bool with_pseudo) const
{
    if (!with_pseudo)
        return _ret_params->mass_var;
    return log_sum(_ret_params->mass_var, 2*(_alpha+_cached_eff_len-2));
}

double Target::log_likelihood(const FragHit& frag, bool with_pseudo) const
{
    double ll = 0;
    
    if (_globs->mismatch_table)
        ll += (_globs->mismatch_table)->log_likelihood(frag);
    
    const PairStatus ps = frag.pair_status();
    {        
        if (with_pseudo)
            ll += log_sum(_ret_params->mass, _alpha+_cached_eff_len + _avg_bias);
        else
            ll += _ret_params->mass;
        
        if (_globs->bias_table)
        {
            if (ps != RIGHT_ONLY)
                ll += _start_bias->at(frag.left);
            if (ps != LEFT_ONLY)
                ll += _end_bias->at(frag.right-1);  
        }
        ll -= (_cached_eff_len + _avg_bias);
    }
    
    if (ps == PAIRED)
        ll += (_globs->fld)->pdf(frag.length());

    assert(!(isnan(ll)||isinf(ll)));
    return ll;
}

double Target::est_effective_length(FLD* fld, bool with_bias) const
{
    if (!fld)
    {
        fld = _globs->fld;
    }
    
    double eff_len = HUGE_VAL;
    
    for(size_t l = fld->min_val(); l <= min(length(), fld->max_val()); l++)
    {
        eff_len = log_sum(eff_len, fld->pdf(l)+log((double)length()-l+1));
    }
    
    
    if (with_bias)
    {
        eff_len += _avg_bias;
    }
    return eff_len;
}

double Target::cached_effective_length(bool with_bias) const
{
    if (with_bias)
        return _cached_eff_len + _avg_bias;
    return _cached_eff_len;
}


void Target::update_target_bias(BiasBoss* bias_table, FLD* fld)
{
    if (!bias_table)
    {
        bias_table = _globs->bias_table;
        fld = _globs->fld;
    }
    else
    {
        _avg_bias = bias_table->get_target_bias(*_start_bias, *_end_bias, *this);
        assert(!isnan(_avg_bias) && !isinf(_avg_bias));
        
    }
    double eff_len = est_effective_length(fld, false);
    {
        _cached_eff_len = eff_len;
    }
}

TargetTable::TargetTable(const string& targ_fasta_file, const TransIndex& targ_index, const TransIndex& targ_lengths, bool prob_seqs, double alpha, const AlphaMap* alpha_map, Globals* globs)
: _globs(globs),
  _targ_map(targ_index.size(), NULL),
  _total_fpb(log(alpha*targ_index.size()))
{
    cout << "Loading target sequences";
    if (globs->bias_table)
        cout << " and measuring bias background";
    cout << "...\n\n";
    
    boost::unordered_set<string> target_names;
    
    double alpha_renorm = 1.0;
    if (alpha_map)
    {
        double alpha_total = 0;
        for(AlphaMap::const_iterator it = alpha_map->begin(); it != alpha_map->end(); ++it)
            alpha_total += it->second;
        alpha_renorm = (alpha * alpha_map->size())/alpha_total;
    }
    
    ifstream infile (targ_fasta_file.c_str());
    string line;
    string seq = "";
    string name = "";
    if (infile.is_open())
    {
        while ( infile.good() )
        {
            getline (infile, line, '\n');
            if (line[0] == '>')
            {
                if (!name.empty())
                {
                    add_targ(name, seq, prob_seqs, (alpha_map) ? alpha_renorm * alpha_map->find(name)->second : alpha, targ_index, targ_lengths);
                }
                name = line.substr(1,line.find(' ')-1);
                if (target_names.count(name))
                {
                    cerr << "ERROR: Target '" << name << "' is duplicated in the input FASTA.  Ensure all target names are unique and re-map before re-running eXpress\n";
                    exit(1);
                }
                if (alpha_map && !alpha_map->count(name))
                {
                    cerr << "ERROR: Target '" << name << "' is was not found in the prior parameter file.\n";
                    exit(1);
                }

                target_names.insert(name);
                seq = "";
            }
            else
            {
                seq += line;
            }
            
        }
        if (!name.empty())
        {
            add_targ(name, seq, prob_seqs, (alpha_map) ? alpha_renorm * alpha_map->find(name)->second : alpha, targ_index, targ_lengths);
        }
        infile.close();
        if (globs->bias_table)
        {
            globs->bias_table->normalize_expectations();
        }
    }
    else 
    {
        cerr << "ERROR: Unable to open MultiFASTA file '" << targ_fasta_file << "'.\n" ; 
        exit(1);
    }
    
    if (size() == 0)
    {
        cerr << "ERROR: No targets found in MultiFASTA file '" << targ_fasta_file << "'.\n" ; 
        exit(1);        
    }
    
    for(TransIndex::const_iterator it = targ_index.begin(); it != targ_index.end(); ++it)
    {
        if (!_targ_map[it->second])
        {
            cerr << "ERROR: Sequence for target '" << it->first << "' not found in MultiFasta file '" << targ_fasta_file << "'.\n";
            exit(1);
        }
    }
}

TargetTable::~TargetTable()
{
    foreach( Target* targ, _targ_map)
    {
        delete targ;
    }
}

void TargetTable::add_targ(const string& name, const string& seq, bool prob_seq, double alpha, const TransIndex& targ_index, const TransIndex& targ_lengths)
{
    TransIndex::const_iterator it = targ_index.find(name);
    if(it == targ_index.end())
    {
        cerr << "Warning: Target '" << name << "' exists in MultiFASTA but not alignment (SAM/BAM) file.\n";
        return;
    }
    
    if (targ_lengths.find(name)->second != seq.length())
    {
        cerr << "ERROR: Target '" << name << "' differs in length between MultiFASTA and alignment (SAM/BAM) files ("<< seq.length() << " vs. " << targ_lengths.find(name)->second << ").\n";
        exit(1);
    }

    Target* targ = new Target(it->second, name, seq, prob_seq, alpha, _globs);
    if (_globs->bias_table)
        (_globs->bias_table)->update_expectations(*targ);
    _targ_map[targ->id()] = targ;
    targ->bundle(_bundle_table.create_bundle(targ));
}

Target* TargetTable::get_targ(TargID id)
{
    return _targ_map[id];
}

Bundle* TargetTable::merge_bundles(Bundle* b1, Bundle* b2)
{
    if (b1 != b2)
    {
        return _bundle_table.merge(b1, b2);
    }
    return b1;
}

size_t TargetTable::num_bundles()
{
    return _bundle_table.size();
}

void TargetTable::round_reset()
{
    foreach(Target* targ, _targ_map)
    {
        targ->round_reset();
    }
}

void project_to_polytope(vector<Target*> bundle_targ, vector<double>& targ_counts, double bundle_counts)
{
    vector<bool> polytope_bound(bundle_targ.size(), false);
    
    while (true)
    {
        double unbound_counts = 0;
        double bound_counts = 0;
        for (size_t i = 0; i < bundle_targ.size(); ++i)
        {
            Target& targ = *bundle_targ[i];
            
            if (targ_counts[i] > targ.tot_counts())
            {
                targ_counts[i] = targ.tot_counts();
                polytope_bound[i] = true;
            }
            else if (targ_counts[i] < targ.uniq_counts())
            {
                targ_counts[i] = targ.uniq_counts();
                polytope_bound[i] = true;
            }
            
            if (polytope_bound[i])
            {
                bound_counts += targ_counts[i];
            }
            else
            {
                unbound_counts += targ_counts[i];
            }
        }
        
        if (approx_eq(unbound_counts + bound_counts, bundle_counts))
            return;
        
        double normalizer = (bundle_counts - bound_counts)/unbound_counts;
        
        bool unbound_exist = false;
        for (size_t i = 0; i < bundle_targ.size(); ++i)
        {    
            if (!polytope_bound[i])
            {
                targ_counts[i] *= normalizer;
                unbound_exist = true;
            }
        }
                
        if (!unbound_exist)
            polytope_bound = vector<bool>(bundle_targ.size(), false);
    }
}

void TargetTable::output_results(string output_dir, size_t tot_counts, bool output_varcov, bool output_edits)
{ 
    FILE * expr_file = fopen((output_dir + "/results.xprs").c_str(), "w");
    ofstream varcov_file;
    ofstream edits_file;
    if (output_varcov)
        varcov_file.open((output_dir + "/varcov.xprs").c_str());    
    
    if (output_edits)
    {
        edits_file.open((output_dir + "/edits.xprs").c_str());
        edits_file << "target_id\tposition\tp_value\tref_nuc\tP(A)\tP(C)\tP(G)\tP(T)\tobs_A\tobs_C\tobs_G\tobs_T\texp_A\texp_C\texp_G\texp_T\n";
    }
    
    fprintf(expr_file, "bundle_id\ttarget_id\tlength\teff_length\ttot_counts\tuniq_counts\test_counts\teff_counts\tambig_distr_alpha\tambig_distr_beta\tfpkm\tfpkm_conf_low\tfpkm_conf_high\tsolvable\n");

    double l_bil = log(1000000000.);
    double l_tot_counts = log((double)tot_counts);
    
    size_t bundle_id = 0;
    foreach (Bundle* bundle, _bundle_table.bundles())
    {
        ++bundle_id;
        
        const vector<Target*>& bundle_targ = bundle->targets();
        
        if (output_varcov)
        {
            varcov_file << ">" << bundle_id << ": ";
            for (size_t i = 0; i < bundle_targ.size(); ++i)
            {
                if (i)
                    varcov_file << ", ";
                varcov_file << bundle_targ[i]->name();
            }
            varcov_file << endl;
        }        
        
        // Calculate total counts for bundle and bundle-level rho
        double l_bundle_mass = HUGE_VAL;
        for (size_t i = 0; i < bundle_targ.size(); ++i)
        {
            l_bundle_mass = log_sum(l_bundle_mass, bundle_targ[i]->mass()); 
        }
        
        if (bundle->counts())
        {
            double l_bundle_counts = log((double)bundle->counts());
            double l_var_renorm = 2*(l_bundle_counts - l_bundle_mass);
            
            vector<double> targ_counts(bundle_targ.size(),0);
            bool requires_projection = false;

            for (size_t i = 0; i < bundle_targ.size(); ++i)
            {
                Target& targ = *bundle_targ[i];
                double l_targ_frac = targ.mass() - l_bundle_mass;
                targ_counts[i] = sexp(l_targ_frac + l_bundle_counts);
                if (targ_counts[i] > (double)targ.tot_counts() || targ_counts[i] < (double)targ.uniq_counts())
                    requires_projection = true;
            }
                        
            if (bundle_targ.size() > 1 && requires_projection)
            {
                project_to_polytope(bundle_targ, targ_counts, bundle->counts());
            }
            
            // Calculate individual counts and rhos
            for (size_t i = 0; i < bundle_targ.size(); ++i)
            {
                Target& targ = *bundle_targ[i];
                double l_eff_len = targ.est_effective_length();

                // Calculate count variance
                double mass_var = targ.mass_var(false);
                double count_alpha = 0;
                double count_beta = 0;
                double count_var = 0;
                
                if (targ.tot_counts() != targ.uniq_counts())
                {
                    
                    double n = targ.tot_counts()-targ.uniq_counts();
                    double m = (targ_counts[i] - targ.uniq_counts())/n;
                    double v = sexp(targ.var_sum() - targ.tot_ambig_mass());                    
                    
                    double a = -m*(m*m - m + v)/v;
                    double b = (m-1)*(m*m - m + v)/v;
                    if (!targ.solvable())
                    {
                        a = 1;
                        b = 1;
                    }
                    
                    if (targ.solvable() && (v == 0 || a < 0 || b < 0))
                        count_var = mass_var;
                    else
                        count_var = n*a*b*(a+b+n)/((a+b)*(a+b)*(a+b+1));
                    
                    count_alpha = a;
                    count_beta = b;
                    assert(!isnan(count_var) && !isinf(count_var));
                }
                
                double fpkm_std_dev = sexp(0.5*(mass_var + l_var_renorm));
                double fpkm_constant = sexp(l_bil - l_eff_len - l_tot_counts);
                double targ_fpkm = targ_counts[i] * fpkm_constant;
                double fpkm_lo = max(0.0, (targ_counts[i] - 2*fpkm_std_dev) * fpkm_constant);
                double fpkm_hi = (targ_counts[i] + 2*fpkm_std_dev) * fpkm_constant;
                
                double eff_len = sexp(l_eff_len);
                double eff_counts = targ_counts[i] / eff_len * targ.length();
                
                fprintf(expr_file, "" SIZE_T_FMT "\t%s\t" SIZE_T_FMT "\t%f\t" SIZE_T_FMT "\t" SIZE_T_FMT "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%c\n", bundle_id, targ.name().c_str(), targ.length(), eff_len, targ.tot_counts(), targ.uniq_counts(), targ_counts[i], eff_counts, count_alpha, count_beta, targ_fpkm, fpkm_lo, fpkm_hi, (targ.solvable())?'T':'F');
            
                if (output_varcov)
                {
                    for (size_t j = 0; j < bundle_targ.size(); ++j)
                    {
                        if (j)
                            varcov_file << "\t";
                        
                        if (i==j)
                            varcov_file << scientific << count_var;
                        else
                            varcov_file << scientific << -sexp(get_covar(targ.id(), bundle_targ[j]->id()) + l_var_renorm);
                           
                    }
                    varcov_file << endl;
                }
                if (output_edits)
                {
                    const Sequence& targ_seq = targ.seq();
                    vector<double> p_vals;
                    targ_seq.calc_p_vals(p_vals);
                    for (size_t i = 0; i < p_vals.size(); ++i)
                    {
                        if (p_vals[i] < 0.05)
                        {
                            edits_file << targ.name() << "\t" << i << "\t" << p_vals[i] << "\t" << NUCS[targ_seq.get_ref(i)];
                            for (size_t nuc=0; nuc < NUM_NUCS; nuc++)
                                edits_file << "\t" << sexp(targ_seq.get_prob(i,nuc));
                            for (size_t nuc=0; nuc < NUM_NUCS; nuc++)
                                edits_file << "\t" << sexp(targ_seq.get_obs(i,nuc));
                            for (size_t nuc=0; nuc < NUM_NUCS; nuc++)
                                edits_file << "\t" << sexp(targ_seq.get_exp(i,nuc));
                            edits_file << endl; 
                        }
                    }
                }
            }
        }
        else
        {
            for (size_t i = 0; i < bundle_targ.size(); ++i)
            {
                Target& targ = *bundle_targ[i];
                fprintf(expr_file, "" SIZE_T_FMT "\t%s\t" SIZE_T_FMT "\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%c\n", bundle_id, targ.name().c_str(), targ.length(), sexp(targ.est_effective_length()), 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 'T');
                
                if (output_varcov)
                {
                    for (size_t j = 0; j < bundle_targ.size(); ++j)
                    {
                        if (j)
                            varcov_file << "\t";
                        varcov_file << scientific << 0.0;
                    }
                    varcov_file << endl;
                }
            }   
        }

    }
    fclose(expr_file);
    if (output_varcov)
        varcov_file.close();
    if (output_edits)
        edits_file.close();
}

double TargetTable::total_fpb() const
{
    boost::unique_lock<boost::mutex>(_fpb_mut);
    return _total_fpb;
}

void TargetTable::update_total_fpb(double incr_amt)
{
    boost::unique_lock<boost::mutex>(_fpb_mut);
    _total_fpb = log_sum(_total_fpb, incr_amt);
}

void TargetTable::threaded_bias_update(boost::mutex* mut)
{
    BiasBoss* bias_table = NULL;
    BiasBoss* bg_table = NULL;
    FLD* fld = NULL;
    
    bool burned_out_before = false;
    
    while(running)
    {
        
        if (bg_table)
            bg_table->normalize_expectations();
        {
            boost::unique_lock<boost::mutex> lock(*mut);
            if(!fld)
                fld = new FLD(*(_globs->fld));
            else
                *fld = *(_globs->fld);
            
            if (_globs->bias_table)
            {
                BiasBoss& glob_bias_table = *(_globs->bias_table);
                if (!bias_table)
                {
                    bias_table = new BiasBoss(glob_bias_table);
                }
                else
                {
                    glob_bias_table.copy_expectations(*bg_table);
                    bg_table->copy_observations(glob_bias_table);                    
                    delete bias_table;
                    bias_table = bg_table;
                }
                bg_table = new BiasBoss(0);
            }    
            cout << "Synchronized parameter tables.\n";
        }

        if (!edit_detect && burned_out && burned_out_before)
            break;
        
        burned_out_before = burned_out;
        
        vector<double> fl_cdf = fld->cdf();
        
        foreach(Target* targ, _targ_map)
        {  
            targ->lock();
            targ->update_target_bias(bias_table, fld);
            if (bg_table)
            {
                bg_table->update_expectations(*targ, targ->rho(), fl_cdf);
            }
            targ->unlock();
            if (!running)
                break;
        }
    }
    
    if (fld)
        delete fld;
    if (bias_table)
        delete bias_table;
    if (bg_table)
        delete bg_table;
}

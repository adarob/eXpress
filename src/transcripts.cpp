//
//  transcripts.cpp
//  express
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "main.h"
#include "transcripts.h"
#include "fld.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cassert>
#include <stdio.h>

using namespace std;

Transcript::Transcript(const TransID id, const std::string& name, const std::string& seq, double alpha, const Globals* globs)
:   _globs(globs),
    _id(id),
    _name(name),
    _seq_f(seq, 0),
    _seq_r(seq, 1),
    _alpha(log(alpha)),
    _ret_params(&_curr_params),
    _uniq_counts(0),
    _tot_counts(0),
    _avg_bias(0),
    _solveable(false)
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

void Transcript::add_mass(double p, double v, double mass) 
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
    
    (_globs->trans_table)->update_total_fpb(mass - _cached_eff_len);
    
}  

void Transcript::round_reset()
{
    _last_params = _curr_params;
    _curr_params = RoundParams();
    _ret_params = &_last_params;
}

double Transcript::rho() const
{
    double eff_len = cached_effective_length(false);
    if (eff_len == HUGE_VAL)
        return HUGE_VAL;
    
    return mass(true) - eff_len - (_globs->trans_table)->total_fpb();
}

double Transcript::mass(bool with_pseudo) const
{
    if (!with_pseudo)
        return _ret_params->mass;
    return log_sum(_ret_params->mass, _alpha+_cached_eff_len);
}

double Transcript::mass_var(bool with_pseudo) const
{
    if (!with_pseudo)
        return _ret_params->mass_var;
    return log_sum(_ret_params->mass_var, 2*(_alpha+_cached_eff_len-2));
}

double Transcript::log_likelihood(const FragHit& frag, bool with_pseudo) const
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

double Transcript::est_effective_length(FLD* fld, bool with_bias) const
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

double Transcript::cached_effective_length(bool with_bias) const
{
    if (with_bias)
        return _cached_eff_len + _avg_bias;
    return _cached_eff_len;
}


void Transcript::update_transcript_bias(BiasBoss* bias_table, FLD* fld)
{
    if (!bias_table)
    {
        bias_table = _globs->bias_table;
        fld = _globs->fld;
    }
    else
    {
        _avg_bias = bias_table->get_transcript_bias(*_start_bias, *_end_bias, *this);
        assert(!isnan(_avg_bias) && !isinf(_avg_bias));
        
    }
    double eff_len = est_effective_length(fld, false);
    {
        _cached_eff_len = eff_len;
    }
}

TranscriptTable::TranscriptTable(const string& trans_fasta_file, const TransIndex& trans_index, const TransIndex& trans_lengths, double alpha, const AlphaMap* alpha_map, Globals* globs)
: _globs(globs),
  _trans_map(trans_index.size(), NULL),
  _total_fpb(log(alpha*trans_index.size()))
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
    
    ifstream infile (trans_fasta_file.c_str());
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
                    add_trans(name, seq, (alpha_map) ? alpha_renorm * alpha_map->find(name)->second : alpha, trans_index, trans_lengths);
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
            add_trans(name, seq, (alpha_map) ? alpha_renorm * alpha_map->find(name)->second : alpha, trans_index, trans_lengths);
        }
        infile.close();
    }
    else 
    {
        cerr << "ERROR: Unable to open MultiFASTA file '" << trans_fasta_file << "'.\n" ; 
        exit(1);
    }
    
    if (size() == 0)
    {
        cerr << "ERROR: No targets found in MultiFASTA file '" << trans_fasta_file << "'.\n" ; 
        exit(1);        
    }
    
    for(TransIndex::const_iterator it = trans_index.begin(); it != trans_index.end(); ++it)
    {
        if (!_trans_map[it->second])
        {
            cerr << "ERROR: Sequence for target '" << it->first << "' not found in MultiFasta file '" << trans_fasta_file << "'.\n";
            exit(1);
        }
    }
}

TranscriptTable::~TranscriptTable()
{
    foreach( Transcript* trans, _trans_map)
    {
        delete trans;
    }
}

void TranscriptTable::add_trans(const string& name, const string& seq, double alpha, const TransIndex& trans_index, const TransIndex& trans_lengths)
{
    TransIndex::const_iterator it = trans_index.find(name);
    if(it == trans_index.end())
    {
        cerr << "Warning: Target '" << name << "' exists in MultiFASTA but not alignment (SAM/BAM) file.\n";
        return;
    }
    
    if (trans_lengths.find(name)->second != seq.length())
    {
        cerr << "ERROR: Target '" << name << "' differs in length between MultiFASTA and alignment (SAM/BAM) files ("<< seq.length() << " vs. " << trans_lengths.find(name)->second << ").\n";
        exit(1);
    }

    Transcript* trans = new Transcript(it->second, name, seq, alpha, _globs);
    if (_globs->bias_table)
        (_globs->bias_table)->update_expectations(*trans);
    _trans_map[trans->id()] = trans;
    trans->bundle(_bundle_table.create_bundle(trans));
}

Transcript* TranscriptTable::get_trans(TransID id)
{
    return _trans_map[id];
}

Bundle* TranscriptTable::merge_bundles(Bundle* b1, Bundle* b2)
{
    if (b1 != b2)
    {
        return _bundle_table.merge(b1, b2);
    }
    return b1;
}

size_t TranscriptTable::num_bundles()
{
    return _bundle_table.size();
}

void TranscriptTable::round_reset()
{
    foreach(Transcript* trans, _trans_map)
    {
        trans->round_reset();
    }
}

void project_to_polytope(vector<Transcript*> bundle_trans, vector<double>& trans_counts, double bundle_counts)
{
    vector<bool> polytope_bound(bundle_trans.size(), false);
    
    while (true)
    {
        double unbound_counts = 0;
        double bound_counts = 0;
        for (size_t i = 0; i < bundle_trans.size(); ++i)
        {
            Transcript& trans = *bundle_trans[i];
            
            if (trans_counts[i] > trans.tot_counts())
            {
                trans_counts[i] = trans.tot_counts();
                polytope_bound[i] = true;
            }
            else if (trans_counts[i] < trans.uniq_counts())
            {
                trans_counts[i] = trans.uniq_counts();
                polytope_bound[i] = true;
            }
            
            if (polytope_bound[i])
            {
                bound_counts += trans_counts[i];
            }
            else
            {
                unbound_counts += trans_counts[i];
            }
        }
        
        if (unbound_counts + bound_counts == bundle_counts)
            return;
        
        double normalizer = (bundle_counts - bound_counts)/unbound_counts;
        bool unbound_exist = false;
        unbound_counts = 0;
        for (size_t i = 0; i < bundle_trans.size(); ++i)
        {    
            if (!polytope_bound[i])
            {
                trans_counts[i] *= normalizer;
                unbound_counts += trans_counts[i];
                unbound_exist = true;
            }
        }
        
        if (unbound_counts + bound_counts - bundle_counts < EPSILON)
            return;
        
        if (!unbound_exist)
            polytope_bound = vector<bool>(bundle_trans.size(), false);
    }
}

void TranscriptTable::output_results(string output_dir, size_t tot_counts, bool output_varcov)
{ 
    FILE * expr_file = fopen((output_dir + "/results.xprs").c_str(), "w");
    ofstream varcov_file;
    if (output_varcov)
        varcov_file.open((output_dir + "/varcov.xprs").c_str());    
    
    fprintf(expr_file, "bundle_id\ttarget_id\tlength\teff_length\ttot_counts\tuniq_counts\test_counts\teff_counts\tambig_distr_alpha\tambig_distr_beta\tfpkm\tfpkm_conf_low\tfpkm_conf_high\tsolveable\n");

    double l_bil = log(1000000000.);
    double l_tot_counts = log((double)tot_counts);
    
    size_t bundle_id = 0;
    foreach (Bundle* bundle, _bundle_table.bundles())
    {
        ++bundle_id;
        
        const vector<Transcript*>& bundle_trans = bundle->transcripts();
        
        if (output_varcov)
        {
            varcov_file << ">" << bundle_id << ": ";
            for (size_t i = 0; i < bundle_trans.size(); ++i)
            {
                if (i)
                    varcov_file << ", ";
                varcov_file << bundle_trans[i]->name();
            }
            varcov_file << endl;
        }        
        
        // Calculate total counts for bundle and bundle-level rho
        double l_bundle_mass = HUGE_VAL;
        for (size_t i = 0; i < bundle_trans.size(); ++i)
        {
            l_bundle_mass = log_sum(l_bundle_mass, bundle_trans[i]->mass()); 
        }
        
        if (bundle->counts())
        {
            double l_bundle_counts = log((double)bundle->counts());
            double l_var_renorm = 2*(l_bundle_counts - l_bundle_mass);
            
            vector<double> trans_counts(bundle_trans.size(),0);
            bool requires_projection = false;

            for (size_t i = 0; i < bundle_trans.size(); ++i)
            {
                Transcript& trans = *bundle_trans[i];
                double l_trans_frac = trans.mass() - l_bundle_mass;
                trans_counts[i] = sexp(l_trans_frac + l_bundle_counts);
                if (!approx_eq(trans_counts[i], (double)trans.tot_counts()))
                    requires_projection = true;
            }
            
            
            if (bundle_trans.size() > 1 && requires_projection)
            {
                project_to_polytope(bundle_trans, trans_counts, bundle->counts());
            }
            
            
            // Calculate individual counts and rhos
            for (size_t i = 0; i < bundle_trans.size(); ++i)
            {
                Transcript& trans = *bundle_trans[i];
                double l_eff_len = trans.est_effective_length();

                // Calculate count variance
                double mass_var = trans.mass_var(false);
                double count_alpha = 0;
                double count_beta = 0;
                double count_var = 0;
                
                if (trans.tot_counts() != trans.uniq_counts())
                {
                    
                    double n = trans.tot_counts()-trans.uniq_counts();
                    double m = (trans_counts[i] - trans.uniq_counts())/n;
                    double v = sexp(trans.var_sum() - trans.tot_ambig_mass());                    
                    
                    double a = -m*(m*m - m + v)/v;
                    double b = (m-1)*(m*m - m + v)/v;
                    if (!trans.solveable())
                    {
                        a = 1;
                        b = 1;
                    }
                    
                    if (trans.solveable() && (v == 0 || a < 0 || b < 0))
                        count_var = mass_var;
                    else
                        count_var = n*a*b*(a+b+n)/((a+b)*(a+b)*(a+b+1));
                    
                    count_alpha = a;
                    count_beta = b;
                    assert(!isnan(count_var) && !isinf(count_var));
                }
                
                double fpkm_std_dev = sexp(0.5*(mass_var + l_var_renorm));
                double fpkm_constant = sexp(l_bil - l_eff_len - l_tot_counts);
                double trans_fpkm = trans_counts[i] * fpkm_constant;
                double fpkm_lo = max(0.0, (trans_counts[i] - 2*fpkm_std_dev) * fpkm_constant);
                double fpkm_hi = (trans_counts[i] + 2*fpkm_std_dev) * fpkm_constant;
                
                double eff_len = sexp(l_eff_len);
                double eff_counts = trans_counts[i] / eff_len * trans.length();
                
                fprintf(expr_file, "" SIZE_T_FMT "\t%s\t" SIZE_T_FMT "\t%f\t" SIZE_T_FMT "\t" SIZE_T_FMT "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%c\n", bundle_id, trans.name().c_str(), trans.length(), eff_len, trans.tot_counts(), trans.uniq_counts(), trans_counts[i], eff_counts, count_alpha, count_beta, trans_fpkm, fpkm_lo, fpkm_hi, (trans.solveable())?'T':'F');
            
                if (output_varcov)
                {
                    for (size_t j = 0; j < bundle_trans.size(); ++j)
                    {
                        if (j)
                            varcov_file << "\t";
                        
                        if (i==j)
                            varcov_file << scientific << count_var;
                        else
                            varcov_file << scientific << -sexp(get_covar(trans.id(), bundle_trans[j]->id()) + l_var_renorm);
                           
                    }
                    varcov_file << endl;
                }
            }
        }
        else
        {
            for (size_t i = 0; i < bundle_trans.size(); ++i)
            {
                Transcript& trans = *bundle_trans[i];
                fprintf(expr_file, "" SIZE_T_FMT "\t%s\t" SIZE_T_FMT "\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%c\n", bundle_id, trans.name().c_str(), trans.length(), sexp(trans.est_effective_length()), 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 'T');
                
                if (output_varcov)
                {
                    for (size_t j = 0; j < bundle_trans.size(); ++j)
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
}

double TranscriptTable::total_fpb() const
{
    boost::unique_lock<boost::mutex>(_fpb_mut);
    return _total_fpb;
}

void TranscriptTable::update_total_fpb(double incr_amt)
{
    boost::unique_lock<boost::mutex>(_fpb_mut);
    _total_fpb = log_sum(_total_fpb, incr_amt);
}

void TranscriptTable::threaded_bias_update(boost::mutex* mut)
{
    BiasBoss* bias_table = NULL;
    BiasBoss* bg_table = NULL;
    FLD* fld = NULL;
    
    bool burned_out_before = false;
    
    while(running)
    {
        
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

        if (burned_out && burned_out_before)
            break;
        
        burned_out_before = burned_out;
        
        vector<double> fl_cdf = fld->cdf();
        
        foreach(Transcript* trans, _trans_map)
        {  
            trans->lock();
            trans->update_transcript_bias(bias_table, fld);
            if (bg_table)
            {
                bg_table->update_expectations(*trans, trans->rho(), fl_cdf);
            }
            trans->unlock();
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

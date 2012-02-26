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
    _avg_bias(0)
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
}

void Transcript::add_mass(double p, double v, double mass) 
{ 
    boost::mutex::scoped_lock lock(_bias_lock);
    _curr_params.mass = log_sum(_curr_params.mass, p+mass);
    _curr_params.samp_var = log_sum(_curr_params.samp_var, 2*mass+p);
    if (p != 0.0 || v != HUGE_VAL)
    {
        _curr_params.binom_var = log_sum(_curr_params.binom_var, 2*mass + p + log(1-sexp(p)));
        _curr_params.tot_unc = log_sum(_curr_params.tot_unc, v+mass);
        _curr_params.ambig_mass = log_sum(_curr_params.ambig_mass, mass+p);
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
    boost::mutex::scoped_lock lock(_bias_lock);
    if (!with_pseudo)
        return _ret_params->mass;
    return log_sum(_ret_params->mass, _alpha+_cached_eff_len);
}

double Transcript::log_likelihood(const FragHit& frag, bool with_pseudo) const
{
    double ll = 0;
    
    if (_globs->mismatch_table)
        ll += (_globs->mismatch_table)->log_likelihood(frag);
    
    const PairStatus ps = frag.pair_status();
    {
        boost::mutex::scoped_lock lock(_bias_lock);
        
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
        boost::mutex::scoped_lock lock(_bias_lock);
        eff_len += _avg_bias;
    }
    return eff_len;
}

double Transcript::cached_effective_length(bool with_bias) const
{
    boost::mutex::scoped_lock lock(_bias_lock);
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
        boost::mutex::scoped_lock lock(_bias_lock);
        _avg_bias = bias_table->get_transcript_bias(*_start_bias, *_end_bias, *this);
        assert(!isnan(_avg_bias) && !isinf(_avg_bias));
        
    }
    double eff_len = est_effective_length(fld, false);
    {
        boost::mutex::scoped_lock lock(_bias_lock);
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
        if (globs->bias_table)
        {
            globs->bias_table->normalize_expectations();
        }
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
    
    fprintf(expr_file, "bundle_id\ttarget_id\tlength\teff_length\ttot_counts\tuniq_counts\tpost_count_mode\tpost_count_var\teff_count_mode\teff_count_var\tfpkm\tfpkm_conf_low\tfpkm_conf_high\n");

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
            l_bundle_mass = log_sum(l_bundle_mass, bundle_trans[i]->mass(false)); 
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
                double l_trans_frac = trans.mass(false) - l_bundle_mass;
                trans_counts[i] = sexp(l_trans_frac + l_bundle_counts);
                if (trans_counts[i] - (double)trans.tot_counts() > EPSILON ||  (double)trans.uniq_counts() - trans_counts[i] > EPSILON)
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
                double count_var = 0;
                
                
                if (trans.tot_counts() != trans.uniq_counts())
                {
                    double binom_var = min(sexp(trans.binom_var() + l_var_renorm), 0.25*trans.tot_counts());
                    
                    long double m = sexp(trans.ambig_mass() - trans.tot_ambig_mass());
                    long double v = sexp(trans.tot_uncertainty() - trans.tot_ambig_mass());
                    //assert (p >=0 && p <= 1);
                    double n = trans.tot_counts()-trans.uniq_counts();
                    
                    long double m2 = m*m;
                    long double m3 = m2*m;
                    long double m4 = m3*m;
                    long double m5 = m4*m;                    
                    long double m6 = m5*m;
                    long double m7 = m6*m;
                    long double m8 = m7*m;
                    long double m9 = m8*m;
                    
                    long double v2 = v*v;
                    long double v3 = v2*v;
                    
                    long double c = -m3 + 2*m2 - 7*v*m - m + 4*v;
                    long double d = 2*m9 - 12*m8 + 42*v*m7 + 30*m7 - 210*v*m6 - 40*m6 + 150*v2*m5 + 429*v*m5 + 30*m5 - 600*v2*m4 - 456*v*m4 - 12*m4 + 2*v3*m3 + 936*v2*m3 + 264*v*m3 + 2*m3 - 6*v3*m2 - 708*v2*m2 - 78*v*m2 + 6*v3*m + 258*v2*m + 9*v*m - 2*v3 - 36*v2;
                    long double e = 2*m3 + 16*v*m2 - 5*m2 - 18*v*m + 4*m + 5*v - 1;
                    long double f = 3*v*e - c*c;
                    long double x = 4*pow(f,3)+d*d;
                    long double y = 4*f*f*f + d*d;
                    if (x < 0)
                        x = 0;
                    if (y < 0)
                        y = 0;
                    long double b = -c/(3*v) + pow(d + sqrt(x), (long double)1.0/3)/(3.77976315*v) - 1.25992105 * f / (3*v*pow(d + sqrt(y), (long double)1.0/3));
                    long double a = (-b*m + 2*m - 1)/(m-1);
                    
                    if (v == 0 || a < 0 || b < 0)
                        count_var = binom_var;
                    else
                        count_var = n*a*b*(a+b+n)/((a+b)*(a+b)*(a+b+1));
                    count_var = binom_var;
                    assert(!isnan(count_var) && !isinf(count_var));
                }
                
                double fpkm_std_dev = sqrt(trans_counts[i] + count_var);
                double fpkm_constant = sexp(l_bil - l_eff_len - l_tot_counts);
                double trans_fpkm = trans_counts[i] * fpkm_constant;
                double fpkm_lo = max(0.0, (trans_counts[i] - 2*fpkm_std_dev) * fpkm_constant);
                double fpkm_hi = (trans_counts[i] + 2*fpkm_std_dev) * fpkm_constant;
                
                double eff_len = sexp(l_eff_len);
                double eff_count_mode = trans_counts[i] / eff_len * trans.length();
                double eff_count_var = count_var * trans.length() * trans.length() / (eff_len*eff_len);
                
                fprintf(expr_file, "" SIZE_T_FMT "\t%s\t" SIZE_T_FMT "\t%f\t" SIZE_T_FMT "\t" SIZE_T_FMT "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", bundle_id, trans.name().c_str(), trans.length(), eff_len, trans.tot_counts(), trans.uniq_counts(), trans_counts[i], count_var, eff_count_mode, eff_count_var, trans_fpkm, fpkm_lo, fpkm_hi);
            
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
                fprintf(expr_file, "" SIZE_T_FMT "\t%s\t" SIZE_T_FMT "\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", bundle_id, trans.name().c_str(), trans.length(), sexp(trans.est_effective_length()), 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                
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

        if (burned_out && burned_out_before)
            break;
        
        burned_out_before = burned_out;
        
        foreach(Transcript* trans, _trans_map)
        {  
            trans->update_transcript_bias(bias_table, fld);
            if (bg_table)
            {
                bg_table->update_expectations(*trans, trans->rho());
            }
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

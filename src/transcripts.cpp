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
#include <assert.h>
#include <stdio.h>

using namespace std;


Transcript::Transcript(const TransID id, const std::string& name, const std::string& seq, double alpha, const Globals* globs)
:   _globs(globs),
    _id(id),
    _name(name),
    _seq(seq),
    _mass_var(HUGE_VAL),
    _est_counts(0),
    _est_counts_var(0),
    _uniq_counts(0),
    _tot_counts(0),
    _start_bias(std::vector<double>(seq.length(),0)),
    _end_bias(std::vector<double>(seq.length(),0)),
    _avg_bias(0)
{   
    _ub_eff_len = est_effective_length();
    _mass = log(_ub_eff_len*alpha); 
}

void Transcript::add_mass(double p, double mass) 
{ 
    assert(sexp(p) != 0);
    _mass = log_sum(_mass, p+mass);
    _tot_counts++;
    if (p != 0.0)
    {
        _mass_var = log_sum(_mass_var, 2*mass + p + log(1-sexp(p)));
    }
}  

void Transcript::add_prob_count(double p)
{
    _est_counts += p;
    _est_counts_var += p*(1-p);
}

double Transcript::log_likelihood(const FragHit& frag) const
{
    boost::mutex::scoped_lock lock(_bias_lock);

    double ll = mass();
    if (_globs->mismatch_table)
        ll += (_globs->mismatch_table)->log_likelihood(frag);
    
    switch(frag.pair_status())
    {
        case PAIRED:
        {
            ll += _start_bias[frag.left] + _end_bias[frag.right-1];
            ll += (_globs->fld)->pdf(frag.length());
            ll -= total_bias_for_length(frag.length());
            break;
        }
        case LEFT_ONLY:
        {
            ll += _start_bias[frag.left];
            ll -= total_bias_for_length(min((int)length()-frag.left, (int)sexp((_globs->fld)->mean())));
            break;
        }
        case RIGHT_ONLY:
        {
            ll += _end_bias[frag.right-1];
            ll -= total_bias_for_length(min(frag.right, (int)sexp((_globs->fld)->mean())));
            break;
        }
    }
    
    return ll;
}

double Transcript::est_effective_length() const
{
    double eff_len = 0.0;
    
    for(size_t l = 1; l <= min(length(), (_globs->fld)->max_val()); l++)
    {
        eff_len += sexp((_globs->fld)->pdf(l))*(length()-l+1);
    }
    
    boost::mutex::scoped_lock lock(_bias_lock);
    eff_len *= sexp(_avg_bias);
    return eff_len;
}

double Transcript::unbiased_effective_length() const
{
    return _ub_eff_len;
}

double Transcript::effective_length() const
{
    double eff_len = 0.0;
    boost::mutex::scoped_lock lock(_bias_lock);
    
    for(size_t l = 1; l <= min(length(), (_globs->fld)->max_val()); l++)
    {
        double len_bias = 0;
        for (size_t i = 0; i < length()-l+1; i++)
        {
            len_bias += sexp(_start_bias[i] + _end_bias[i+l]);
        }
        eff_len += sexp((_globs->fld)->pdf(l))*len_bias;
    }
    
    return eff_len;
}

double Transcript::total_bias_for_length(size_t l) const
{
    assert(l <= (_globs->fld)->max_val());
    return _avg_bias + log((double)(length() - l + 1));
}

void Transcript::update_transcript_bias()
{
    if (!_globs->bias_table)
        return;
    boost::mutex::scoped_lock lock(_bias_lock);
    _avg_bias = (_globs->bias_table)->get_transcript_bias(_start_bias, _end_bias, *this);
}

TranscriptTable::TranscriptTable(const string& trans_fasta_file, const TransIndex& trans_index, double alpha, const Globals* globs)
: _globs(globs),
  _trans_map(trans_index.size(), NULL),
  _alpha(alpha)
{
    cout << "Loading target sequences from " << trans_fasta_file << "...\n\n";
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
                    add_trans(name, seq, trans_index);
                }
                name = line.substr(1,line.find(' ')-1);
                seq = "";
            }
            else
            {
                seq += line;
            }
            
        }
        if (!name.empty())
        {
            add_trans(name, seq, trans_index);
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

void TranscriptTable::add_trans(const string& name, const string& seq, const TransIndex& trans_index)
{
    TransIndex::const_iterator it = trans_index.find(name);
    if(it == trans_index.end())
    {
        cerr << "Warning: Target '" << name << "' exists in MultiFASTA but not alignment (SAM/BAM) file.\n";
        return;
    }
    Transcript* trans = new Transcript(it->second, name, seq, _alpha, _globs);
    if (_globs->bias_table)
        (_globs->bias_table)->update_expectations(*trans);
    _trans_map[trans->id()] = trans;
    trans->bundle(_bundle_table.create_bundle(trans));
}

Transcript* TranscriptTable::get_trans(TransID id)
{
    return _trans_map[id];
}

void TranscriptTable::update_covar(TransID trans1, TransID trans2, double covar)
{
    size_t pair_id = size()*min(trans1, trans2)+max(trans1, trans2);
    if (_covar_map.count(pair_id))
    {
        _covar_map[pair_id] = log_sum(covar, _covar_map[pair_id]);
    }
    else
    {
        _covar_map[pair_id] = covar;
    }
}

double TranscriptTable::get_covar(TransID trans1, TransID trans2)
{
    size_t pair_id = size()*min(trans1, trans2)+max(trans1, trans2);
    if (_covar_map.count(pair_id))
    {
        return _covar_map[pair_id];
    }
    else
    {
        return HUGE_VAL;
    }
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


void TranscriptTable::threaded_bias_update()
{
    while(running)
    {
        //size_t count = 0;
        foreach(Transcript* trans, _trans_map)
        {  
            trans->update_transcript_bias();
            //if (++count%10000 == 0)
            //    cout << "*"<<count<<"\n";
            if (!running)
                break;
        }
        //if (size() > 1000)
        //    cout << "* Completed bias update of all " << size() << " targets.\n";
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
        
        if (unbound_counts + bound_counts - bundle_counts < 0.00001)
            return;
        
        if (!unbound_exist)
            polytope_bound = vector<bool>(bundle_trans.size(), false);
    }
}

void TranscriptTable::output_results(string output_dir, size_t tot_counts, bool output_varcov, bool multi_iteration)
{
    FILE * expr_file = fopen((output_dir + "/results.xprs").c_str(), "w");
    ofstream varcov_file;
    if (output_varcov)
        varcov_file.open((output_dir + "/varcov.xprs").c_str());    
    
    fprintf(expr_file, "bundle_id\ttarget_id\tlength\teff_length\ttot_counts\tuniq_counts\tpost_count_mean\tpost_count_var\teff_count_mean\teff_count_var\tfpkm\tfpkm_conf_low\tfpkm_conf_high\n");

    double l_bil = log(1000000000.);
    double l_tot_counts = log((double)tot_counts);

    // Bundle transcripts based on partition
    typedef boost::unordered_map<TransID, vector<Transcript*> > BundleMap;
    BundleMap b_map;
    
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
                trans_counts[i] = (multi_iteration) ? trans.est_counts():sexp(l_trans_frac + l_bundle_counts);
            
                if (trans_counts[i] > (double)trans.tot_counts() || trans_counts[i] < (double)trans.uniq_counts())
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
                double eff_len = trans.est_effective_length();

                double count_var = (multi_iteration) ? trans.est_counts_var():min(sexp(trans.mass_var() + l_var_renorm), 0.25*trans.tot_counts());
                double eff_count_norm = (double)trans.length()/eff_len;
                
                double fpkm_std_dev = sqrt(trans_counts[i] + count_var);
                double fpkm_constant = sexp(l_bil - log(eff_len) - l_tot_counts);
                double trans_fpkm = trans_counts[i] * fpkm_constant;
                double fpkm_lo = max(0.0, (trans_counts[i] - 2*fpkm_std_dev) * fpkm_constant);
                double fpkm_hi = (trans_counts[i] + 2*fpkm_std_dev) * fpkm_constant;
                
                fprintf(expr_file, "%zu\t%s\t%zu\t%f\t%zu\t%zu\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", bundle_id, trans.name().c_str(), trans.length(), eff_len, trans.tot_counts(), trans.uniq_counts(), trans_counts[i], count_var, trans_counts[i]*eff_count_norm, count_var*pow(eff_count_norm,2.0),trans_fpkm, fpkm_lo, fpkm_hi);
            
                if (output_varcov)
                {
                    for (size_t j = 0; j < bundle_trans.size(); ++j)
                    {
                        if (j)
                            varcov_file << "\t";
                        
                        if (i==j)
                            varcov_file << scientific << count_var;
                        else if (multi_iteration)
                            varcov_file << scientific << -get_covar(trans.id(), bundle_trans[j]->id());
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
                fprintf(expr_file, "%zu\t%s\t%zu\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", bundle_id, trans.name().c_str(), trans.length(), trans.est_effective_length(), 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                
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

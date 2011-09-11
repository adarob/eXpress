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


Transcript::Transcript(const std::string& name, const std::string& seq, double alpha, const FLD* fld, const BiasBoss* bias_table, const MismatchTable* mismatch_table)
:_name(name),
_id(hash_trans_name(name)),
_seq(seq),
_var(HUGE_VAL),
_uniq_counts(0),
_tot_counts(0),
_bundle_counts(0),
_start_bias(std::vector<double>(seq.length(),0)),
_end_bias(std::vector<double>(seq.length(),0)),
_avg_bias(0),
_fld(fld),
_bias_table(bias_table),
_mismatch_table(mismatch_table)
{   _ub_eff_len = est_effective_length();
    _mass = log(_ub_eff_len*alpha); }

void Transcript::add_mass(double p, double mass) 
{ 
    assert(sexp(p) != 0);
    _mass = log_sum(_mass, p+mass);
    _tot_counts++;
    if (p == 0.0)
        // probability of 1
    {
        _uniq_counts++;
    }
    else
    {
        _var = log_sum(_var, 2*mass + p + log(1-sexp(p)));
    }
}  

double Transcript::log_likelihood(const FragMap& frag) const
{
    boost::mutex::scoped_lock lock(_bias_lock);

    double ll = mass();
    if (_mismatch_table)
        ll += _mismatch_table->log_likelihood(frag);
    
    switch(frag.pair_status())
    {
        case PAIRED:
        {
            ll += _start_bias[frag.left] + _end_bias[frag.right-1];
            ll += _fld->pdf(frag.length());
            ll -= total_bias_for_length(frag.length());
            break;
        }
        case LEFT_ONLY:
        {
            ll += _start_bias[frag.left];
            ll -= total_bias_for_length(min((int)length()-frag.left, (int)sexp(_fld->mean())));
            break;
        }
        case RIGHT_ONLY:
        {
            ll += _end_bias[frag.right-1];
            ll -= total_bias_for_length(min(frag.right, (int)sexp(_fld->mean())));
            break;
        }
    }
    
    return ll;
}

double Transcript::est_effective_length() const
{
    double eff_len = 0.0;
    
    for(size_t l = 1; l <= min(length(), _fld->max_val()); l++)
    {
        eff_len += sexp(_fld->pdf(l))*(length()-l+1);
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
    
    for(size_t l = 1; l <= min(length(), _fld->max_val()); l++)
    {
        double len_bias = 0;
        for (size_t i = 0; i < length()-l+1; i++)
        {
            len_bias += sexp(_start_bias[i] + _end_bias[i+l]);
        }
        eff_len += sexp(_fld->pdf(l))*len_bias;
    }
    
    boost::mutex::scoped_lock lock(_bias_lock);
    return eff_len;
}

double Transcript::total_bias_for_length(size_t l) const
{
    assert(l <= _fld->max_val());
    return _avg_bias + log((double)(length() - l + 1));
}

void Transcript::update_transcript_bias()
{
    if (!_bias_table)
        return;
    boost::mutex::scoped_lock lock(_bias_lock);
    _avg_bias = _bias_table->get_transcript_bias(_start_bias, _end_bias, *this);
}

TranscriptTable::TranscriptTable(const string& trans_fasta_file, double alpha, const FLD* fld, BiasBoss* bias_table, const MismatchTable* mismatch_table)
: _rank(Rank(_rank_map)),
  _parent(Parent(_parent_map)),
  _bundles(_rank, _parent),
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
                    Transcript* trans = new Transcript(name, seq, alpha, fld, bias_table, mismatch_table);
                    add_trans(trans);
                    if (bias_table)
                        bias_table->update_expectations(*trans);
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
            Transcript* trans = new Transcript(name, seq, alpha, fld, bias_table, mismatch_table);
            add_trans(trans);
            if (bias_table)
                bias_table->update_expectations(*trans);
        }
        infile.close();
    }
    else 
    {
        cerr << "Unable to open MultiFASTA file '" << trans_fasta_file << "'.\n" ; 
        exit(1);
    }
    
    if (size() == 0)
    {
        cerr << "No targets found in MultiFASTA file '" << trans_fasta_file << "'.\n" ; 
        exit(1);        
    }
    
}

TranscriptTable::~TranscriptTable()
{
    for( TransMap::iterator it = _trans_map.begin(); it != _trans_map.end(); ++it)
    {
        delete it->second;
    }
}

void TranscriptTable::add_trans(Transcript* trans)
{
    Transcript* ret = get_trans(trans->id());
    if (ret)
    {
        if (trans->name() == ret->name())
        {
            cerr << "ERROR: Repeated target ID in multi-fasta input (" << ret->name() << ").\n";
        }
        else
        {
            cerr << "ERROR: Hash collision (" << ret->name() << ").\n";
        }
        exit(1);
    }
    
    _trans_map.insert(make_pair(trans->id(), trans));
    _bundles.make_set(trans->id());
}

Transcript* TranscriptTable::get_trans(TransID id)
{
    TransMap::iterator it = _trans_map.find(id);
    if(it != _trans_map.end())
        return it->second;
    return NULL;
}

void TranscriptTable::update_covar(TransID trans1, TransID trans2, double covar)
{
    size_t pair_id = hash_trans_pair(trans1, trans2);
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
    size_t pair_id = hash_trans_pair(trans1, trans2);
    if (_covar_map.count(pair_id))
    {
        return _covar_map[pair_id];
    }
    else
    {
        return HUGE_VAL;
    }
}


TransID TranscriptTable::get_trans_rep(TransID trans)
{
    return _bundles.find_set(trans);
}


TransID TranscriptTable::merge_bundles(TransID rep1, TransID rep2)
{
    assert(rep1 != rep2);
    _bundles.link(rep1, rep2);
    
    TransID new_rep = _bundles.find_set(rep1);
    return new_rep;
}

size_t TranscriptTable::num_bundles()
{
    size_t count = 0;
    for(TransMap::const_iterator it = _trans_map.begin(); it != _trans_map.end(); ++it)
    {
        if (_bundles.find_set(it->first) == it->first)
            count ++;
    }
    return count;
}

void TranscriptTable::output_bundles(string output_dir)
{
    ofstream bundle_file((output_dir + "/bundles.tab").c_str());
    
    typedef boost::unordered_map<TransID, string> BundleMap;
    BundleMap b_map;
    for( TransMap::iterator it = _trans_map.begin(); it != _trans_map.end(); ++it)
    {
        TransID rep = _bundles.find_set(it->first);
        b_map[rep] += '\t' + (it->second)->name();
    }
    
    for (BundleMap::iterator it = b_map.begin(); it != b_map.end(); ++it)
    {
        bundle_file << it->second << '\n';
    }
    
    bundle_file.close();
}

void TranscriptTable::threaded_bias_update()
{
    while(running)
    {
        //size_t count = 0;
        for( TransMap::iterator it = _trans_map.begin(); it != _trans_map.end(); ++it)
        {  
            Transcript& trans = *(it->second);
            trans.update_transcript_bias();
            //if (++count%10000 == 0)
            //    cout << "*"<<count<<"\n";
            if (!running)
                break;
        }
        //if (size() > 1000)
        //    cout << "* Completed bias update of all " << size() << " targets.\n";
    }
}

void TranscriptTable::output_header(ofstream& runexpr_file)
{
  runexpr_file << "Read_Num\t";
    for( TransMap::iterator it = _trans_map.begin(); it != _trans_map.end(); ++it)
    {
        Transcript& trans = *(it->second);
        runexpr_file << trans.name() << '\t'; 
    }   
    runexpr_file << '\n';
}

void TranscriptTable::output_current(ofstream& runexpr_file)
{    
    double sum = HUGE_VAL;
    for( TransMap::iterator it = _trans_map.begin(); it != _trans_map.end(); ++it)
    {
        Transcript& trans = *(it->second);
        double log_fpkm = trans.mass() - log(trans.est_effective_length());
        sum = log_sum(sum, log_fpkm);
    }   
    
    for( TransMap::iterator it = _trans_map.begin(); it != _trans_map.end(); ++it)
    {

        Transcript& trans = *(it->second);
        double log_fpkm = trans.mass() - log(trans.est_effective_length());
        runexpr_file << sexp(log_fpkm - sum) << '\t';
    }   
    runexpr_file << '\n';
}

void project_to_polytope(vector<Transcript*>& bundle_trans, vector<double>& trans_counts, vector<bool>& polytope_bound, double unbound_counts, double req_unbound_counts)
{
    double new_bound_counts = 0;
    double new_unbound_counts = 0;
    bool new_fixed = false;
    for (size_t i = 0; i < bundle_trans.size(); ++i)
    {
        if (polytope_bound[i])
            continue;
        
        Transcript& trans = *bundle_trans[i];
        trans_counts[i] *= req_unbound_counts/unbound_counts;
        
        if (trans_counts[i] > trans.tot_counts())
        {
            trans_counts[i] = trans.tot_counts();
            new_bound_counts += trans_counts[i];
            polytope_bound[i] = true;
            new_fixed = true;
        }
        else if (trans_counts[i] < trans.uniq_counts())
        {
            trans_counts[i] = trans.uniq_counts();
            new_bound_counts += trans_counts[i];
            polytope_bound[i] = true;
            new_fixed = true;
        }
        else
        {
            trans_counts[i] *= req_unbound_counts/unbound_counts;
            new_unbound_counts += trans_counts[i];
        }
    }
    
    if (new_fixed)
    {
        project_to_polytope(bundle_trans, trans_counts, polytope_bound, new_unbound_counts, req_unbound_counts - new_bound_counts);
    }
    return;
}

void TranscriptTable::output_results(string output_dir, size_t tot_counts, bool output_varcov)
{
    FILE * expr_file = fopen((output_dir + "/results.xprs").c_str(), "w");
    ofstream varcov_file;
    if (output_varcov)
        varcov_file.open((output_dir + "/varcov.xprs").c_str());    
    
    fprintf(expr_file, "bundle_id\ttarget_id\tlength\teff_length\ttot_counts\tuniq_counts\test_counts\test_counts_var\tfpkm\tfpkm_conf_low\tfpkm_conf_high\n");
    double l_bil = log(1000000000.);
    double l_tot_counts = log((double)tot_counts);

    // Bundle transcripts based on partition
    typedef boost::unordered_map<TransID, vector<Transcript*> > BundleMap;
    BundleMap b_map;
    for( TransMap::iterator it = _trans_map.begin(); it != _trans_map.end(); ++it)
    {
        TransID rep = _bundles.find_set(it->first);
        b_map[rep].push_back(it->second);
    }
    
    size_t bundle_id = 0;
    for (BundleMap::iterator it = b_map.begin(); it != b_map.end(); ++it)
    {
        ++bundle_id;
        vector<Transcript*> bundle_trans = it->second;
        
        if (output_varcov)
            varcov_file << ">" << bundle_id << ": ";
        
        // Calculate total counts for bundle and bundle-level rho
        size_t bundle_counts = 0;
        double l_bundle_mass = HUGE_VAL;
        for (size_t i = 0; i < bundle_trans.size(); ++i)
        {
            bundle_counts += bundle_trans[i]->bundle_counts();
            l_bundle_mass = log_sum(l_bundle_mass, bundle_trans[i]->mass()); 
            if (output_varcov)
            {
                if (i)
                    varcov_file << ", ";
                varcov_file << bundle_trans[i]->name();
            }
        }
        if (output_varcov)
            varcov_file << endl;
        
        if (bundle_counts)
        {
            double l_bundle_counts = log((double)bundle_counts);
            double l_var_renorm = 2*(l_bundle_counts - l_bundle_mass);
            
            vector<double> trans_counts(bundle_trans.size(),0);
            bool requires_projection = false;

            for (size_t i = 0; i < bundle_trans.size(); ++i)
            {
                Transcript& trans = *bundle_trans[i];
                double l_trans_frac = trans.mass() - l_bundle_mass;
                trans_counts[i] = sexp(l_trans_frac + l_bundle_counts);
            
                if (trans_counts[i] > trans.tot_counts() || trans_counts[i] < trans.uniq_counts())
                    requires_projection = true;
            }
            
            if (requires_projection)
            {
                vector<bool> polytope_bound(bundle_trans.size(),false);
                project_to_polytope(bundle_trans, trans_counts, polytope_bound, bundle_counts, bundle_counts);
            }
            
            
            // Calculate individual counts and rhos
            for (size_t i = 0; i < bundle_trans.size(); ++i)
            {
                Transcript& trans = *bundle_trans[i];
                double l_trans_var = trans.var() + l_var_renorm;
                double std_dev = sexp(0.5*l_trans_var);
                double eff_len = trans.est_effective_length();
                double fpkm_constant = sexp(l_bil - log(eff_len) - l_tot_counts);
                double trans_fpkm = trans_counts[i] * fpkm_constant;
                double fpkm_lo = max(0.0, (trans_counts[i] - 2*std_dev) * fpkm_constant);
                double fpkm_hi = (trans_counts[i] + 2*std_dev) * fpkm_constant;
                fprintf(expr_file, "%zu\t%s\t%zu\t%f\t%zu\t%zu\t%f\t%f\t%f\t%f\t%f\n", bundle_id, trans.name().c_str(), trans.length(), eff_len, trans.tot_counts(), trans.uniq_counts(), trans_counts[i], sexp(l_trans_var), trans_fpkm, fpkm_lo, fpkm_hi);
            
                if (output_varcov)
                {
                    for (size_t j = 0; j < bundle_trans.size(); ++j)
                    {
                        if (j)
                            varcov_file << " ";
                        if (i==j)
                            varcov_file << scientific << sexp(l_trans_var);
                        else
                        {
                            varcov_file << scientific << -sexp(get_covar(trans.id(), bundle_trans[j]->id()) + l_var_renorm);
                        }   
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
                fprintf(expr_file, "%zu\t%s\t%zu\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n", bundle_id, trans.name().c_str(), trans.length(), trans.effective_length(), 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0);
                
                if (output_varcov)
                {
                    for (size_t j = 0; j < bundle_trans.size(); ++j)
                    {
                        if (j)
                            varcov_file << " ";
                        varcov_file << scientific << 0.0;
                    }
                    varcov_file << endl;
                }
            }   
        }

    }
}

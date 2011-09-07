//
//  transcripts.cpp
//  express
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "main.h"
#include "transcripts.h"
#include "bundles.h"
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
{   
    _ub_eff_len = est_effective_length();
    _mass = log(_ub_eff_len*alpha); 
}

void Transcript::add_mass(double p, double mass) 
{ 
    assert(sexp(p) != 0);
    _mass = log_sum(_mass, p+mass);
    _tot_counts++;
    if (p == 0.0)
        // probabilitiy of 1
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

TranscriptTable::TranscriptTable(const string& trans_fasta_file, double alpha, const FLD* fld, FragMassTable* fmt, BiasBoss* bias_table, const MismatchTable* mismatch_table)
: _alpha(alpha),
  _bundle_table(fmt)
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
    trans->bundle(_bundle_table.create_bundle(trans));
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


void TranscriptTable::output_results(string output_dir, size_t tot_counts, bool output_varcov)
{
    FILE * expr_file = fopen((output_dir + "/results.xprs").c_str(), "w");
    ofstream varcov_file;
    if (output_varcov)
        varcov_file.open((output_dir + "/varcov.xprs").c_str());    
    
    fprintf(expr_file, "bundle_id\ttarget_id\tlength\teff_length\ttot_counts\tuniq_counts\test_counts\test_counts_var\tfpkm\tfpkm_conf_low\tfpkm_conf_high\n");
    double l_bil = log(1000000000.);
    double l_tot_counts = log((double)tot_counts);
    
    size_t bundle_id = 0;
    foreach (Bundle* bundle, _bundle_table.bundles())
    {
        ++bundle_id;
        
        vector<Transcript*>& bundle_trans = bundle->transcripts();
        
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
        
        if (bundle->counts())
        {
            double l_bundle_mass = bundle->mass();
            double l_bundle_counts = log((double)bundle->counts());
            double l_var_renorm = 2*(l_bundle_counts - l_bundle_mass);
            // Calculate individual counts and rhos
            for (size_t i = 0; i < bundle_trans.size(); ++i)
            {
                Transcript& trans = *bundle_trans[i];
                double l_trans_frac = trans.mass() - l_bundle_mass;
                double trans_counts = sexp(l_trans_frac + l_bundle_counts);
                double l_trans_var = trans.var() + l_var_renorm;
                double std_dev = sexp(0.5*l_trans_var);
                double eff_len = trans.est_effective_length();
                double fpkm_constant = sexp(l_bil - log(eff_len) - l_tot_counts);
                double trans_fpkm = trans_counts * fpkm_constant;
                double fpkm_lo = max(0.0, (trans_counts - 2*std_dev) * fpkm_constant);
                double fpkm_hi = (trans_counts + 2*std_dev) * fpkm_constant;
                fprintf(expr_file, "%zu\t%s\t%zu\t%f\t%zu\t%zu\t%f\t%f\t%f\t%f\t%f\n", bundle_id, trans.name().c_str(), trans.length(), eff_len, trans.tot_counts(), trans.uniq_counts(), trans_counts, sexp(l_trans_var), trans_fpkm, fpkm_lo, fpkm_hi);
            
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

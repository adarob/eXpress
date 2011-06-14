//
//  transcripts.cpp
//  expressionline
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "transcripts.h"
#include "fld.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>


using namespace std;


Transcript::Transcript(const std::string& name, const std::string& seq, double alpha, const FLD* fld, const BiasBoss* bias_table, const MismatchTable* mismatch_table)
:_name(name),
_id(hash_trans_name(name.c_str())),
_seq(seq),
_len(seq.length()),
_counts(seq.length()*alpha),
_start_bias(std::vector<double>(seq.length(),1.0)),
_end_bias(std::vector<double>(seq.length(),1.0)),
_avg_bias(1),
_fld(fld),
_bias_table(bias_table),
_mismatch_table(mismatch_table)
{ }

double Transcript::log_likelihood(const FragMap& frag) const
{
    boost::mutex::scoped_lock lock(_bias_lock);

    double ll = log(frag_count());
    ll += _mismatch_table->log_likelihood(frag, *this);
    
    switch(frag.pair_status())
    {
        case PAIRED:
        {
            double len_prob = _fld->pdf(frag.length());
            if (len_prob == 0) return 0.0;
            ll += log(_start_bias[frag.left] * _end_bias[frag.right-1]);
            ll += log(len_prob);
            ll -= log(total_bias_for_length(frag.length()));
            break;
        }
        case LEFT_ONLY:
        {
            ll += log(_start_bias[frag.left]);
            ll -= log(total_bias_for_length(min((int)length()-frag.left, (int)_fld->mean())));
            break;
        }
        case RIGHT_ONLY:
        {
            ll += log(_end_bias[frag.right-1]);
            ll -= log(total_bias_for_length(min(frag.right, (int)_fld->mean())));
            break;
        }
    }
    
    assert(!isnan(ll) && !isinf(ll));
    return ll;
}

double Transcript::effective_length() const
{
    double eff_len = 0.0;
    
    for(size_t l = 1; l <= min(length(), _fld->max_val()); l++)
    {
        eff_len += _fld->pdf(l)*(length()-l+1);
    }
    
    boost::mutex::scoped_lock lock(_bias_lock);
    eff_len *= _avg_bias;
    return eff_len;
}

double Transcript::total_bias_for_length(size_t l) const
{
    assert(l <= _fld->max_val());
    return _avg_bias * (length() - l + 1);
}

void Transcript::update_transcript_bias()
{
    if (!_bias_table)
        return;
    boost::mutex::scoped_lock lock(_bias_lock);
    _avg_bias = _bias_table->get_transcript_bias(_start_bias, _end_bias, *this);
}

TranscriptTable::TranscriptTable(const string& trans_fasta_file, double alpha, const FLD* fld, BiasBoss* bias_table, const MismatchTable* mismatch_table)
{
    _alpha = alpha;
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
        cerr << "No transcripts found in MultiFASTA file '" << trans_fasta_file << "'.\n" ; 
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
            cerr << "ERROR: Repeated transcript ID in multi-fasta input (" << ret->name() << ").\n";
        }
        else
        {
            cerr << "ERROR: Hash collision (" << ret->name() << ").\n";
        }
        exit(1);
    }
    
    _trans_map.insert(make_pair(trans->id(), trans));
}

Transcript* TranscriptTable::get_trans(TransID id)
{
    TransMap::iterator it = _trans_map.find(id);
    if(it != _trans_map.end())
        return it->second;
    return NULL;
}

void TranscriptTable::threaded_bias_update()
{
    while(true)
    {
        size_t count = 0;
        for( TransMap::iterator it = _trans_map.begin(); it != _trans_map.end(); ++it)
        {  
            Transcript& trans = *(it->second);
            trans.update_transcript_bias();
            if (++count%1000 == 0)
                cout << "*"<<count<<"\n";

        }
        cout << "*\n";
    }
}

void TranscriptTable::output_expression(string output_dir)
{
    ofstream expr_file((output_dir + "/transcripts.expr").c_str());
    for( TransMap::iterator it = _trans_map.begin(); it != _trans_map.end(); ++it)
    {
        Transcript& trans = *(it->second);
        double M = trans.fld()->num_obs();
        trans.update_transcript_bias();
        double counts = trans.frag_count() - _alpha * trans.length();
        double fpkm = (counts/trans.effective_length())*(1000000000/M);
        expr_file << trans.name() << "\t" << fpkm << "\n";
    }   
    expr_file.close();
}

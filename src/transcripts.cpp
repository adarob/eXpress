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

Transcript::Transcript(std::string& id, const std::string& seq, double alpha, const FLD* fld, const BiasBoss* bias_table, const MismatchTable* mismatch_table)
:_id(id),
_seq(seq),
_len(seq.length()),
_counts(seq.length()*alpha),
_start_bias(std::vector<double>(seq.length(),1.0)),
_end_bias(std::vector<double>(seq.length(),1.0)),
_fld(fld),
_bias_table(bias_table),
_mismatch_table(mismatch_table)
{
    _tot_bias_for_len = std::vector<double>(_fld->max_val()+1);
    for (size_t l = 1; l <= _fld->max_val(); ++l)
    {
        _tot_bias_for_len[l] = seq.length()-l+1;
    }
}

double Transcript::likelihood(const FragMap& frag) const
{
    double len_prob = _fld->pdf(frag.length());
    if (len_prob == 0) return 0.0;

    double l = frag_count();
    l *= _start_bias[frag.left] * _end_bias[frag.right-1];
    l *= len_prob;
    l *= _mismatch_table->likelihood(frag, *this);
    l /= total_bias_for_length(frag.length());
    assert(!isnan(l) && !isinf(l));
    return l;
}

double Transcript::log_likelihood(const FragMap& frag) const
{
    double len_prob = _fld->pdf(frag.length());
    if (len_prob == 0) return 0.0;
    
    double ll = log(frag_count());
    ll += log(_start_bias[frag.left] * _end_bias[frag.right-1]);
    ll += log(len_prob);
    ll += _mismatch_table->log_likelihood(frag, *this);
    ll -= log(total_bias_for_length(frag.length()));
    assert(!isnan(ll) && !isinf(ll));
    return ll;
}

double Transcript::effective_length() const
{
    double eff_len = 0.0;
    
    for(size_t l = 1; l <= min(length(), _fld->max_val()); l++)
    {
        eff_len += _fld->pdf(l)*total_bias_for_length(l);
    }
    return eff_len;
}

double Transcript::total_bias_for_length(size_t l) const
{
    //boost::mutex::scoped_lock lock(_bias_lock);
    assert(l <= _fld->max_val());
    return _tot_bias_for_len[l];
}

void Transcript::update_transcript_bias()
{
    if (!_bias_table)
        return;
    boost::mutex::scoped_lock lock(_bias_lock);
    _bias_table->get_transcript_bias(_start_bias, _end_bias, *this);
    for (size_t l = 1; l <= min(length(), _fld->max_val()); ++l)
    {
        double len_bias = 0;
        for(size_t i = 0; i <= length()-l+1; ++i)
        {
            len_bias += _start_bias[i]*_end_bias[i+l-1];
        }
        _tot_bias_for_len[l] = len_bias;
    }
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
                if (name != "")
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
        if (name != "")
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

// This is FNV-1, see http://en.wikipedia.org/wiki/Fowler_Noll_Vo_hash
inline uint64_t TranscriptTable::hash_id(const string& id) const
{
    const char * __s = id.c_str();
    uint64_t hash = 0xcbf29ce484222325ull;
    for ( ; *__s; ++__s)
    {
        hash *= 1099511628211ull;
        hash ^= *__s;
    }
    return hash;
}

void TranscriptTable::add_trans(Transcript* trans)
{
    Transcript* ret = get_trans(trans->id());
    if (ret)
    {
        if (trans->id() == ret->id())
        {
            cerr << "ERROR: Repeated transcript ID in multi-fasta input (" << ret->id() << ").\n";
        }
        else
        {
            cerr << "ERROR: Hash collision (" << ret->id() << ").\n";
        }
        exit(1);
    }
    _trans_map.insert(make_pair(hash_id(trans->id()), trans));
}

Transcript* TranscriptTable::get_trans(const string& name)
{
    TransMap::iterator it = _trans_map.find(hash_id(name));
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
        double M = trans.fld()->total();
        trans.update_transcript_bias();
        double counts = trans.frag_count() - _alpha * trans.length();
        double fpkm = (counts/trans.effective_length())*(1000000000/M);
        expr_file << trans.id() << "\t" << fpkm << "\n";
    }   
    expr_file.close();
}

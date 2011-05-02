//
//  mismatchmodel.cpp
//  expressionline2
//
//  Created by Adam Roberts on 4/23/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "mismatchmodel.h"
#include "transcripts.h"
#include "fragments.h"
#include <iostream>
#include <fstream>

using namespace std;

const size_t MAX_READ_LEN = 200;

inline size_t MismatchTable::ctoi(const char c) const
{
    switch(c)
    {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return 4;
    }
}

inline size_t MismatchTable::ctoi_r(const char c) const
{
    switch(c)
    {
        case 'A':
        case 'a':
            return 3;
        case 'C':
        case 'c':
            return 2;
        case 'G':
        case 'g':
            return 1;
        case 'T':
        case 't':
            return 0;
        default:
            return 4;
    }
}

MismatchTable::MismatchTable(double alpha)
{
    _markov_model = vector<FrequencyMatrix>(MAX_READ_LEN, FrequencyMatrix(16, 4, alpha));
}

double MismatchTable::likelihood(const FragMap& f, const Transcript& t) const
{
    const string& t_seq = t.seq();
    double p = 1.0;
    
    size_t cur;
    size_t ref;
    size_t prev = 0;
    size_t index;
    
    for (size_t i = 0; i < f.seq_l.length(); i++)
    {
        cur = ctoi(f.seq_l[i]);
        ref = ctoi(t_seq[f.left+i]);
        if (prev != 4 && cur != 4 && ref != 4)
        {
            index = (prev << 2) + ref;
            p *= _markov_model[i](index, cur);
        }
        prev = ref;
    }
    
    size_t r_len = f.seq_r.length();
    prev = 0;
    for (size_t i = 1; i < r_len; i++)
    {
        ref = ctoi_r(t_seq[f.right-i-1]);
        cur = ctoi_r(f.seq_r[r_len-i-1]);
        if (prev != 4 && cur != 4 && ref != 4)
        {
            index = (prev << 2) + ref;
            p *= _markov_model[i](index, cur);
        }
        prev = ref;
    }
    return p;
}

double MismatchTable::log_likelihood(const FragMap& f, const Transcript& t) const
{
    const string& t_seq = t.seq();
    double ll = 0;
    
    size_t cur;
    size_t ref;
    size_t prev = 0;
    size_t index;
    
    for (size_t i = 0; i < f.seq_l.length(); i++)
    {
        cur = ctoi(f.seq_l[i]);
        ref = ctoi(t_seq[f.left+i]);
        if (prev != 4 && cur != 4 && ref != 4)
        {
            index = (prev << 2) + ref;
            ll += log(_markov_model[i](index, cur));
        }
        prev = ref;
    }
    
    size_t r_len = f.seq_r.length();
    prev = 0;
    for (size_t i = 1; i < r_len; i++)
    {
        ref = ctoi_r(t_seq[f.right-i-1]);
        cur = ctoi_r(f.seq_r[r_len-i-1]);
        if (prev != 4 && cur != 4 && ref != 4)
        {
            index = (prev << 2) + ref;
            ll += log(_markov_model[i](index, cur));
        }
        prev = ref;
    }
    return ll;
}


void MismatchTable::update(const FragMap& f, const Transcript& t, double mass)
{
    const string& t_seq = t.seq();
    
    size_t cur;
    size_t ref;
    size_t prev = 0;
    size_t index;
    
    for (size_t i = 0; i < f.seq_l.length(); i++)
    {
        cur = ctoi(f.seq_l[i]);
        ref = ctoi(t_seq[f.left+i]);
        if (prev != 4 && cur != 4 && ref != 4)
        {
            index = (prev << 2) + ref;
            _markov_model[i].increment(index, cur, mass);
        }
        prev = ref;
    }
    
    size_t r_len = f.seq_r.length();
    prev = 0;
    for (size_t i = 0; i < r_len; i++)
    {
        ref = ctoi_r(t_seq[f.right-i-1]);
        cur = ctoi_r(f.seq_r[r_len-i-1]);
        if (prev != 4 && cur != 4 && ref != 4)
        {
            index = (prev << 2) + ref;
            _markov_model[i].increment(index, cur, mass);
        }
        prev = ref;
    }
}

void MismatchTable::output(string path)
{
    string filename = path + "/mismatch_probs.tab";
    ofstream outfile(filename.c_str());
 
    for(size_t i = 0; i < 16; i++)
    {
        for(size_t j = 0; j < 4; j++)
        {
            for (size_t k = 0; k < MAX_READ_LEN; k++)
            {
                outfile << _markov_model[k](i,j)<<"\t";
            }
            outfile<<"\n";
        }
    }
}
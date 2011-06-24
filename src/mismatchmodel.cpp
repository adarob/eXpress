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

MismatchTable::MismatchTable(long double alpha)
{
    _first_read_mm = vector<FrequencyMatrix>(MAX_READ_LEN, FrequencyMatrix(16, 4, alpha));
    _second_read_mm = vector<FrequencyMatrix>(MAX_READ_LEN, FrequencyMatrix(16, 4, alpha));
}

long double MismatchTable::likelihood(const FragMap& f, const Transcript& t) const
{
    const string& t_seq = t.seq();
    long double p = 1.0;
    
    const vector<FrequencyMatrix>& left_mm = (f.left_first) ? _first_read_mm : _second_read_mm;
    const vector<FrequencyMatrix>& right_mm = (!f.left_first) ? _first_read_mm : _second_read_mm;
    
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
            p *= left_mm[i](index, cur);
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
            p *= right_mm[i](index, cur);
        }
        prev = ref;
    }
    return p;
}

long double MismatchTable::log_likelihood(const FragMap& f, const Transcript& t) const
{
    const string& t_seq = t.seq();
    long double ll = 0;
    
    const vector<FrequencyMatrix>& left_mm = (f.left_first) ? _first_read_mm : _second_read_mm;
    const vector<FrequencyMatrix>& right_mm = (!f.left_first) ? _first_read_mm : _second_read_mm;
    
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
            ll += log(left_mm[i](index, cur));
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
            ll += log(right_mm[i](index, cur));
        }
        prev = ref;
    }
    assert(!isnan(ll) && !isinf(ll));
    return ll;
}


void MismatchTable::update(const FragMap& f, const Transcript& t, long double mass)
{
    const string& t_seq = t.seq();
    
    vector<FrequencyMatrix>& left_mm = (f.left_first) ? _first_read_mm : _second_read_mm;
    vector<FrequencyMatrix>& right_mm = (!f.left_first) ? _first_read_mm : _second_read_mm;
    
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
            left_mm[i].increment(index, cur, mass);
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
            right_mm[i].increment(index, cur, mass);
        }
        prev = ref;
    }
}

void MismatchTable::output(string path)
{
    string filename = path + "/mismatch_probs.tab";
    ofstream outfile(filename.c_str());
    
    outfile<<"First Read\n";
    for(size_t k = 0; k < MAX_READ_LEN; k++)
    {
        for(size_t i = 0; i < 16; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                outfile << _first_read_mm[k](i,j)<<"\t";
            }
        }
        outfile<<"\n";
    }

    
    outfile<<"Second Read\n";
    for(size_t k = 0; k < MAX_READ_LEN; k++)
    {
        for(size_t i = 0; i < 16; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                outfile << _second_read_mm[k](i,j)<<"\t";
            }
        }
        outfile<<"\n";
    }
}
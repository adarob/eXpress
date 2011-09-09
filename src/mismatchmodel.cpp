//
//  mismatchmodel.cpp
//  express
//
//  Created by Adam Roberts on 4/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "main.h"
#include "mismatchmodel.h"
#include "transcripts.h"
#include "fragments.h"
#include <iostream>
#include <fstream>

using namespace std;

MismatchTable::MismatchTable(double alpha)
: _first_read_mm(MAX_READ_LEN, FrequencyMatrix(16, 4, alpha)),
  _second_read_mm(MAX_READ_LEN, FrequencyMatrix(16, 4, alpha)),
  _active(false),
  _max_len(0)
{}

double MismatchTable::log_likelihood(const FragMap& f) const
{
    if (!_active)
        return 0;
    
    const string& t_seq = f.mapped_trans->seq();
    double ll = 0;
    
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
            ll += left_mm[i](index, cur);
            assert(!isnan(ll));
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
            ll += right_mm[i](index, cur);
            assert(!isnan(ll));
        }
        prev = ref;
    }
    return ll;
}


void MismatchTable::update(const FragMap& f, double mass)
{
    const string& t_seq = f.mapped_trans->seq();
    
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
    
    _max_len = max(_max_len, max(f.seq_l.length(), f.seq_r.length())); 
}

string MismatchTable::to_string() const
{
    string s = "";
    char buff[50];

    for (size_t ref = 0; ref < NUM_NUCS; ref++)
    {
        vector<double> mass(4, HUGE_VAL);
        double tot = HUGE_VAL;
        for (size_t prev = 0; prev < NUM_NUCS; prev++)
        {
            size_t col_i = (prev << 2) + ref;
            for (size_t obs = 0; obs < NUM_NUCS; obs++)
            {
                size_t arr_i = (col_i << 2) + obs;
                for (size_t k = 0; k < MAX_READ_LEN; k++)
                {
                    double val =  _first_read_mm[k].arr(arr_i);
                    mass[obs] = log_sum(mass[obs], val);
                    tot = log_sum(tot, val);
                }
            }
        }
        for (size_t obs = 0; obs < NUM_NUCS; obs++)
        {
            sprintf(buff, "%e ", sexp(mass[obs]-tot));
            s += buff;
        }
    }
    
    s.erase(s.length()-1,1);
    return s;
}

void MismatchTable::append_output(ofstream& outfile) const
{
    string col_header =  "\t";
    for(size_t i = 0; i < 64; i++)
    {
        col_header += NUCS[i>>4];
        col_header += NUCS[i>>2 & 3]; 
        col_header += "->*";
        col_header += NUCS[i & 3];
        col_header += '\t';
    }
    col_header[col_header.length()-1] = '\n';
    
    outfile<<">First Read Mismatch\n" << col_header;
    for(size_t k = 0; k < _max_len; k++)
    {
        outfile << k << ":\t";
        for(size_t i = 0; i < 16; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                outfile << scientific << sexp(_first_read_mm[k](i,j))<<"\t";
            }
        }
        outfile<<endl;
    }

    
    outfile<<">Second Read Mismatch\n" << col_header;
    for(size_t k = 0; k < _max_len; k++)
    {
        outfile << k << ":\t";
        for(size_t i = 0; i < 16; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                outfile << scientific << sexp(_second_read_mm[k](i,j))<<"\t";
            }
        }
        outfile<<endl;
    }
}
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
{
    _first_read_mm = vector<FrequencyMatrix>(MAX_READ_LEN, FrequencyMatrix(16, 4, alpha));
    _second_read_mm = vector<FrequencyMatrix>(MAX_READ_LEN, FrequencyMatrix(16, 4, alpha));
}

double MismatchTable::log_likelihood(const FragMap& f) const
{
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
                outfile << scientific << sexp(_first_read_mm[k](i,j))<<"\t";
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
                outfile << scientific << sexp(_second_read_mm[k](i,j))<<"\t";
            }
        }
        outfile<<"\n";
    }
}
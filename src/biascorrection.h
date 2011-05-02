#ifndef BIASCORRECTION_H
#define BIASCORRECTION_H

/*
 *  biascorrection.h
 *  cufflinks
 *
 *  Created by Adam Roberts on 4/5/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 *
 */

#include <vector>
#include <string>
#include <boost/thread.hpp>
#include "frequencymatrix.h"

class Transcript;
class FragMap;

class SeqWeightTable
{
    FrequencyMatrix _observed;
    FrequencyMatrix _expected;
    mutable boost::mutex _lock;
    size_t ctoi(char c) const;
public:
    SeqWeightTable(){}
    SeqWeightTable(size_t window_size, double alpha);
    
    void increment_expected(char c); 
    void increment_observed(std::string& seq, double normalized_mass);
    double get_weight(const std::string& seq, size_t i) const;
};


class BiasBoss
{
    SeqWeightTable _5_seq_bias;
    SeqWeightTable _3_seq_bias;
    
public:
    BiasBoss(){};
    BiasBoss(double alpha);
    
    void update_expectations(const Transcript& trans);
    void update_observed(const FragMap& frag, const Transcript& trans, double mass);
    void get_transcript_bias(std::vector<double>& start_bias, std::vector<double>& end_bias, const Transcript& trans) const;
};


#endif

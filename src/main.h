#ifndef MAIN_H
#define MAIN_H
//
//  main.h
//  express
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "config.h"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <boost/foreach.hpp>

//#include <boost/random/mersenne_twister.hpp>

#define foreach BOOST_FOREACH

class TranscriptTable;
class BiasBoss;
class MismatchTable;
class FLD;

/**
 * a struct for holding pointers to the global parameter tables (bias_table, mismatch_table, fld)
 */
struct Globals
{
    FLD* fld;
    MismatchTable* mismatch_table;
    BiasBoss* bias_table;
    TranscriptTable* trans_table;
};

/**
 * a global pseudo-random number generator
 */ 
//extern boost::mt19937 random_gen;

/**
 * a global bool that is true when processing is still occuring
 * this is primarily used to notify the bias update thread to stop
 */
extern bool running;

/**
 * a global size_t specifying the maximum read length supported
 */
const size_t MAX_READ_LEN = 200;

/**
 * a global size_t specifying the number of possible nucleotides
 */
const size_t NUM_NUCS = 4;

/**
 * a global specifying the nucleotide ordering
 */
const char NUCS[] = {'A','C','G','T'};

const double EPSILON = 0.0001;

/**
 * global function to calculate the log of the sum of 2 logged values efficiently
 * log(0) == HUGE_VAL
 * @param x a double for the first logged value in the sum
 * @param y a double for the second logged value in the sum
 * @return a double for the log of the sum of the exp(x)+exp(y)
 */
inline double log_sum(double x, double y)
{
    if (fabs(x) == HUGE_VAL)
    {
        return y;
    }
    if (fabs(y) == HUGE_VAL)
    {
        return x;
    }
    
    if (y>x)
        std::swap(x,y);
    
    
    double sum = x+log(1+exp(y-x));
    return sum;
}

/**
 * global function to determine if a logged value is 0 in non-log space
 * @param x a double for the logged value to be tested
 * @return True if exp(x)==0 (without rounding)
 */
inline double islzero(double x)
{
    return (fabs(x) == HUGE_VAL);
}

/**
 * global function to exponentiate a logged value
 * @param x a double for the logged value to be exponentiated
 * @return exp(x) or 0 if x = log(0)
 */
inline double sexp(double x)
{
    if (islzero(x))
        return 0.0;
    return exp(x);
}

/**
 * global function to encode a nucleotide character to a size_t value
 * @param c the nucleotide character to be encoded
 * @return a size_t value encoding the nucleotide
 */
inline size_t ctoi(const char c)
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
            return 0;
    }
}

/**
 * global function to convert a nucleotide character to a size_t value for its complement
 * @param c the nucleotide character to be complemented and encoded
 * @return a size_t value encoding the complement of the given nucleotide
 */
inline size_t ctoi_r(const char c)
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
            return 3;
    }
}

inline char complement(char c)
{
    switch(c)
    {
        case 'a':
        case 'A':
            return 'T';
        case 'c':
        case 'C':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
            return 'A';
        default:
            return 'N';
    }
}


#endif
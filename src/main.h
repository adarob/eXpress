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

#define foreach BOOST_FOREACH

class MapParser;
class TargetTable;
class BiasBoss;
class MismatchTable;
class FLD;

/**
 * a global bool that is true when processing is still occuring
 * this is primarily used to notify the bias update thread to stop running
 */
extern bool running;

/**
 * a global bool that is true when the auxilary params are finished burning in
 * this is primarily used to notify the bias update thread to stop updating 
 * certain parameters
 */
extern bool burned_out;
extern bool edit_detect;

enum Direction { BOTH, FR, RF };
extern Direction direction;

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

inline bool approx_eq(double a, double b, double eps=EPSILON)
{
    return abs(a-b) < eps;
}

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

#endif

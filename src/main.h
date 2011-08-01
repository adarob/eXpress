#ifndef MAIN_H
#define MAIN_H
//
//  main.h
//  expressionline
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 UC Berkeley. All rights reserved.
//

#include <cmath>
#include <cassert>
#include <algorithm>


extern bool running;
const size_t MAX_READ_LEN = 200;
const int NUM_NUCS = 4;

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

inline double sexp(double x)
{
    if (fabs(x) == HUGE_VAL)
        return 0.0;
    return exp(x);
}

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
            return 4;
    }
}


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
            return 4;
    }
}



#endif
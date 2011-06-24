#ifndef MAIN_H
#define MAIN_H
//
//  main.h
//  expressionline
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 UC Berkeley. All rights reserved.
//

#include <math.h>

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
    return x+log(1+exp(y-x));
}

inline double sexp(double x)
{
    if (fabs(x) == HUGE_VAL)
        return 0.0;
    return exp(x);
}

#endif
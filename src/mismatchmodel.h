#ifndef MISMATCH_H
#define MISMATCH_H
//
//  mismatchmodel.h
//  expressionline2
//
//  Created by Adam Roberts on 4/23/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "frequencymatrix.h"

class FragMap;
class Transcript;

class MismatchTable
{
    std::vector<FrequencyMatrix> _first_read_mm;
    std::vector<FrequencyMatrix> _second_read_mm;

public:
    MismatchTable(double alpha);
    double likelihood(const FragMap& f, const Transcript& t) const;
    double log_likelihood(const FragMap& f, const Transcript& t) const;
    void update(const FragMap& f, const Transcript& t, double mass);
    std::string to_string() const;
    void output(std::string filename);
};

#endif
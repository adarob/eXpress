//
//  bundles.cpp
//  express
//
//  Created by Adam Roberts on 9/6/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "transcripts.h"
#include "bundles.h"
#include "main.h"

using namespace std;

Bundle::Bundle(Transcript* trans)
: _mass(trans->mass()),
  _counts(trans->tot_counts())
{ _transcripts.push_back(trans); }

void Bundle::renormalize_transcripts(double new_bundle_mass)
{
    double l_mass_renorm = new_bundle_mass - mass();
    double l_var_renorm = 2*(new_bundle_mass - mass());
    
    foreach(Transcript* trans, transcripts())
    {
        trans->mass(trans->mass() + l_mass_renorm);
        trans->var(trans->var() + l_var_renorm);
    }
}

void Bundle::add_mass(double mass)
{
    _mass = log_sum(_mass, mass);
}

void Bundle::incr_counts(size_t incr_amt)
{
    _counts += incr_amt;
}

BundleTable::~BundleTable()
{
    foreach(Bundle* b, _bundles)
    {
        delete b;
    }
}

Bundle* BundleTable::create_bundle(Transcript* trans)
{
    Bundle* b = new Bundle(trans);
    _bundles.insert(b);
    return b;
}

Bundle* BundleTable::merge(Bundle* b1, Bundle* b2)
{
    if (b1==b2)
        return b1;
    
    if (b1->size() < b2->size())
        swap(b1, b2);
    
    b1->renormalize_transcripts(b1->mass()+b2->mass());
    b2->renormalize_transcripts(b1->mass()+b2->mass());
    
    foreach(Transcript* trans, b2->transcripts())
    {
        trans->bundle(b1);
        b1->transcripts().push_back(trans);
    }
    
    b1->incr_counts(b2->counts());
    b1->add_mass(b2->mass());
    
    _bundles.erase(b2);
    delete b2;  
    
    return b1;
}
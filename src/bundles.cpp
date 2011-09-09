//
//  bundles.cpp
//  express
//
//  Created by Adam Roberts on 9/6/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "fragmasstable.h"
#include "transcripts.h"
#include "bundles.h"
#include "main.h"

using namespace std;

Bundle::Bundle(Transcript* trans, FragMassTable* fmt)
: _mass(trans->mass()),
  _counts(trans->tot_counts()),
  _pseudo_mass(trans->mass()),
  _n(0),
  _next_frag_mass(0.0),
  _fmt(fmt)    
{ _transcripts.push_back(trans); }

double Bundle::next_frag_mass()
{
    double frag_mass = _next_frag_mass;
    _next_frag_mass = _fmt->next_frag_mass(++_n, frag_mass);
    return frag_mass;
}

void Bundle::renormalize_transcripts(double mass_norm_const)
{
    double var_norm_const = 2*(mass_norm_const);
    
    foreach(Transcript* trans, transcripts())
    {
        trans->mass(trans->mass() + mass_norm_const);
        trans->var(trans->var() + var_norm_const);
    }
    
    //Still need to do covars!
}

void Bundle::add_mass(double mass)
{
    _mass = log_sum(_mass, mass);
}

void Bundle::add_pseudo_mass(double mass)
{
    _pseudo_mass = log_sum(_pseudo_mass, mass);
}

void Bundle::incr_counts(size_t incr_amt)
{
    _counts += incr_amt;
}

double Bundle::counts_plus_pseudo() const
{
    return _counts + sexp(2*_pseudo_mass - _mass);
}


BundleTable::BundleTable(FragMassTable* fmt)
: _fmt(fmt) {}

BundleTable::~BundleTable()
{
    foreach(Bundle* b, _bundles)
    {
        delete b;
    }
}

Bundle* BundleTable::create_bundle(Transcript* trans)
{
    Bundle* b = new Bundle(trans, _fmt);
    _bundles.insert(b);
    return b;
}

Bundle* BundleTable::merge(Bundle* b1, Bundle* b2)
{
    if (b1==b2)
        return b1;
    
    if (b1->size() < b2->size())
        swap(b1, b2);

    double new_bundle_counts = b1->counts() + b2->counts();
    if (new_bundle_counts)
    {
        double new_bundle_mass;
        double new_frag_mass;
        size_t new_n = _fmt->nearest_stored_mass(new_bundle_counts, new_frag_mass, new_bundle_mass);
        double norm_const = new_bundle_mass - log(b1->counts_plus_pseudo() + b2->counts_plus_pseudo());
        b1->renormalize_transcripts(log(b1->counts_plus_pseudo()) + norm_const - b1->mass());
        b2->renormalize_transcripts(log(b2->counts_plus_pseudo()) + norm_const - b2->mass());
        b1->incr_counts(b2->counts());
        b1->add_pseudo_mass(b2->pseudo_mass());
        b1->mass(new_bundle_mass);
        b1->n(new_n);
        b1->next_frag_mass(new_frag_mass);
    }
    else
    {
        b1->add_pseudo_mass(b2->pseudo_mass());
        b1->add_mass(b2->mass());
    }
    
    foreach(Transcript* trans, b2->transcripts())
    {
        trans->bundle(b1);
        b1->transcripts().push_back(trans);
    }
    
    _bundles.erase(b2);
    delete b2;  
    
    return b1;
}
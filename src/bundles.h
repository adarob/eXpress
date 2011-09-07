//
//  bundles.h
//  express
//
//  Created by Adam Roberts on 9/6/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#ifndef express_bundles_h
#define express_bundles_h

#include <vector>
#include <boost/unordered_set.hpp>


class Transcript;
class FragMassTable;

class Bundle
{
    std::vector<Transcript*> _transcripts;     
    double _mass;
    size_t _counts;
    size_t _n;
    size_t _next_frag_mass;
    FragMassTable* _fmt;
    
public:
    Bundle(Transcript* trans, FragMassTable* fmt);
    
    double next_frag_mass();
    void next_frag_mass(double m) { _next_frag_mass = m; }
    void n(size_t n) { _n = n; }
    
    void renormalize_transcripts(double norm_const);
    void add_mass(double mass);
    void incr_counts(size_t incr_amt=1);
    
    size_t size() const { return _transcripts.size(); }
    std::vector<Transcript*>& transcripts() { return _transcripts; }
    double mass() const { return _mass; }
    void mass(double m) { _mass = m; }
    size_t counts() const { return _counts; }
};


typedef boost::unordered_set<Bundle*> BundleSet;
class BundleTable
{
    BundleSet _bundles;
    FragMassTable* _fmt;

public:
    BundleTable(FragMassTable* fmt);
    ~BundleTable();
    const BundleSet bundles() const { return _bundles; }
    size_t size() const { return _bundles.size(); }
    Bundle* create_bundle(Transcript* trans);
    Bundle* merge(Bundle* b1, Bundle* b2);
};


#endif

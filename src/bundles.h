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

class Bundle
{
    std::vector<Transcript*> _transcripts;     
    double _mass;
    size_t _counts;

public:
    Bundle(Transcript* trans);
    void renormalize_transcripts(double new_bundle_mass);
    void add_mass(double mass);
    void incr_counts(size_t incr_amt=1);

    size_t size() const { return _transcripts.size(); }
    std::vector<Transcript*>& transcripts() { return _transcripts; }
    double mass() const { return _mass; }
    size_t counts() const { return _counts; }
};


typedef boost::unordered_set<Bundle*> BundleSet;
class BundleTable
{
    BundleSet _bundles;

public:
    ~BundleTable();
    const BundleSet bundles() const { return _bundles; }
    size_t size() const { return _bundles.size(); }
    Bundle* create_bundle(Transcript* trans);
    Bundle* merge(Bundle* b1, Bundle* b2);
};


#endif

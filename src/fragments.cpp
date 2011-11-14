//
//  fragments.cpp
//  express
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "fragments.h"
#include "main.h"
#include <string.h>
#include <stdlib.h>
//#include <boost/random/discrete_distribution.hpp>

using namespace std;

Fragment::~Fragment()
{
    for (size_t i = 0; i < num_hits(); i++)
    {
        delete _frag_hits[i];
    }
    
    for (size_t i = 0; i < _open_mates.size(); i++)
    {
        delete _open_mates[i];
    }
}

bool Fragment::add_map_end(FragHit* h)
{
    if (_name.empty())
    {
        _name = h->name;
    }
    else if (_name != h->name)
    {
        return false;
    }
    
    if (h->mate_l >= 0)
    {
        add_open_mate(h);
    }
    else
    // single-end fragment
    {
        _frag_hits.push_back(h);
    }
    
    return true;
}

void Fragment::add_open_mate(FragHit* new_p)
{
    bool found = false;
    
    FragHit& nm = *new_p;
    for( vector<FragHit*>::iterator it = _open_mates.begin(); it != _open_mates.end(); ++it)
    {
        FragHit& om = **it;
        if (nm.trans_id == om.trans_id && nm.mate_l == om.left && om.mate_l == nm.left)
        {
            if (nm.left < om.left)
            {
                assert(nm.seq_r.empty());
                nm.right = om.right;
                nm.seq_r = om.seq_r;
                nm.sam_r = om.sam_r;
                nm.bam_r = om.bam_r;
            }
            else
            {
                assert(nm.seq_l.empty());
                nm.left = om.left;
                nm.seq_l = om.seq_l;
                nm.sam_l = om.sam_l;
                nm.bam_l = om.bam_l;
            }
            
            assert(nm.left_first == om.left_first);
            found = true;
            delete *it;
            _frag_hits.push_back(&nm);
            _open_mates.erase(it);
            break;
        }
    }
    
    if (!found)
        _open_mates.push_back(&nm);
}

const FragHit* Fragment::sample_hit() const
{
    vector<double> probs(_frag_hits.size());
    probs[0] = _frag_hits[0]->probability;
    for (size_t i=1; i < _frag_hits.size(); ++i)
    {
        probs[i] = probs[i-1] + _frag_hits[i]->probability;
    }

//   boost::random::discrete_distribution<> dist(probs);
//   return _frag_hits[dist(random_gen)];
    double r = rand()/double(RAND_MAX)*probs.back();
    size_t i = lower_bound(probs.begin(), probs.end(), r) - probs.begin();
    return _frag_hits[i];
}


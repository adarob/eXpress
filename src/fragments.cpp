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

using namespace std;

Fragment::Fragment(Library* lib) : _lib(lib) {}

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
        if (nm.targ_id == om.targ_id && (size_t)nm.mate_l == om.left && (size_t)om.mate_l == nm.left)
        {
	    if (nm.left < om.left || (nm.left == om.left && nm.seq_r.empty()))
            {
                assert(nm.seq_r.empty());
                nm.right = om.right;
                nm.seq_r = om.seq_r;
                nm.sam_r = om.sam_r;
                nm.bam_r = om.bam_r;
                nm.inserts_r = om.inserts_r;
                nm.deletes_r = om.deletes_r;
            }
            else
            {
                assert(nm.seq_l.empty());
                nm.left = om.left;
                nm.seq_l = om.seq_l;
                nm.sam_l = om.sam_l;
                nm.bam_l = om.bam_l;
                nm.inserts_l = om.inserts_l;
                nm.deletes_l = om.deletes_l;
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

    double r = rand()/double(RAND_MAX)*probs.back();
    size_t i = lower_bound(probs.begin(), probs.end(), r) - probs.begin();
    return _frag_hits[i];
}

bool fraghit_compare(FragHit* h1, FragHit* h2)
{
    return h1->targ_id < h2->targ_id;
}

void Fragment::sort_hits()
{
    sort(_frag_hits.begin(), _frag_hits.end(), fraghit_compare);
}

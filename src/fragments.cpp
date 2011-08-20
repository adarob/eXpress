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
using namespace std;

Fragment::~Fragment()
{
    for (size_t i = 0; i < num_maps(); i++)
    {
        delete _frag_maps[i];
    }
    
    for (size_t i = 0; i < _open_mates.size(); i++)
    {
        delete _open_mates[i];
    }
}

bool Fragment::add_map_end(FragMap* m)
{
    if (_name.empty())
    {
        _name = m->name;
    }
    else if (_name != m->name)
    {
        return false;
    }
    
    if (m->mate_l >= 0)
    {
        add_open_mate(m);
    }
    else
    // single-end fragment
    {
        _frag_maps.push_back(m);
    }
    
    return true;
}

void Fragment::add_open_mate(FragMap* new_p)
{
    bool found = false;
    
    FragMap& nm = *new_p;
    for( vector<FragMap*>::iterator it = _open_mates.begin(); it != _open_mates.end(); ++it)
    {
        FragMap& om = **it;
        if (nm.trans_id == om.trans_id && nm.mate_l == om.left && om.mate_l == nm.left)
        {
            if (nm.left < om.left)
            {
                assert(nm.seq_r.empty());
                nm.right = om.right;
                nm.seq_r = om.seq_r;
            }
            else
            {
                assert(nm.seq_l.empty());
                nm.left = om.left;
                nm.seq_l = om.seq_l;
            }
            
            assert(nm.left_first == om.left_first);
            found = true;
            delete *it;
            _frag_maps.push_back(&nm);
            _open_mates.erase(it);
            break;
        }
    }
    
    if (!found)
        _open_mates.push_back(&nm);
}


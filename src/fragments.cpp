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

Fragment::~Fragment() {
  for (size_t i = 0; i < num_hits(); i++) {
    delete _frag_hits[i];
  }

  for (size_t i = 0; i < _open_mates.size(); i++) {
    delete _open_mates[i];
  }
}

bool Fragment::add_map_end(ReadHit* r)
{
  if (_name.empty()) {
    _name = r->name;
  } else if (_name != r->name) {
    return false;
  }

  if (r->mate_l >= 0) {
    add_open_mate(r);
  } else {  // single-end fragment
    _frag_hits.push_back(new FragHit(r));
  }

  return true;
}

void Fragment::add_open_mate(ReadHit* nm) {
  bool found = false;

  for(vector<ReadHit*>::iterator it = _open_mates.begin();
      it != _open_mates.end(); ++it) {
    ReadHit* om = *it;
    if (nm->targ_id == om->targ_id &&
        (size_t)nm->mate_l == om->left &&
        (size_t)om->mate_l == nm->left &&
        nm->first != om->first &&
        nm->reversed != om->reversed) {
      FragHit* h = NULL;
	    if (nm->left < om->left || (nm->left == om->left && om->reversed)) {
        h = new FragHit(nm, om);
      } else {
        h = new FragHit(om, nm);
      }

      found = true;
      _frag_hits.push_back(h);
      _open_mates.erase(it);
      break;
    }
  }

  if (!found) {
    _open_mates.push_back(nm);
  }
}

const FragHit* Fragment::sample_hit() const {
  vector<double> probs(_frag_hits.size());
  probs[0] = sexp(_frag_hits[0]->params()->posterior);
  for (size_t i=1; i < _frag_hits.size(); ++i) {
    probs[i] = probs[i-1] + sexp(_frag_hits[i]->params()->posterior);
  }

  double r = rand()/double(RAND_MAX)*probs.back();
  size_t i = lower_bound(probs.begin(), probs.end(), r) - probs.begin();
  return _frag_hits[i];
}

bool fraghit_compare(FragHit* h1, FragHit* h2) {
  return h1->target_id() < h2->target_id();
}

void Fragment::sort_hits() {
  sort(_frag_hits.begin(), _frag_hits.end(), fraghit_compare);
}

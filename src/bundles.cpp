//
//  bundles.cpp
//  express
//
//  Created by Adam Roberts on 9/6/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "bundles.h"

#include "main.h"
#include "targets.h"

using namespace std;

void CovarTable::increment(TargID targ1, TargID targ2, double incr_amt) {
  size_t pair_id = size()*min(targ1, targ2)+max(targ1, targ2);
  if (_covar_map.count(pair_id)) {
    _covar_map[pair_id] = log_add(_covar_map[pair_id], incr_amt);
  } else {
    _covar_map[pair_id] = incr_amt;
  }
}

double CovarTable::get(TargID targ1, TargID targ2) {
  size_t pair_id = size()*min(targ1, targ2)+max(targ1, targ2);
  if (_covar_map.count(pair_id)) {
    return _covar_map[pair_id];
  } else {
    return LOG_0;
  }
}

Bundle::Bundle(Target* targ)
    : _counts(targ->tot_counts()),
      _mass(targ->mass(true)) {
  _targets.push_back(targ);
}

void Bundle::incr_counts(size_t incr_amt) {
  _counts += incr_amt;
}

void Bundle::incr_mass(double incr_amt) {
  _mass = log_add(_mass, incr_amt);
}

BundleTable::~BundleTable() {
  foreach(Bundle* bundle, _bundles) {
    delete bundle;
  }
}

Bundle* BundleTable::create_bundle(Target* targ) {
  Bundle* b = new Bundle(targ);
  _bundles.insert(b);
  return b;
}

Bundle* BundleTable::merge(Bundle* b1, Bundle* b2) {
  if (b1==b2) {
    return b1;
  }

  if (b1->size() < b2->size()) {
    swap(b1, b2);
  }

  foreach(Target* targ, b2->_targets) {
    targ->bundle(b1);
    b1->_targets.push_back(targ);
  }

  b1->incr_counts(b2->counts());
  b1->incr_mass(b2->mass());
  _bundles.erase(b2);
  delete b2;

  return b1;
}

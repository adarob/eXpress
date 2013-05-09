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
      _mass(targ->mass(true)),
      _merged_into(NULL) {
  _targets.push_back(targ);
}

size_t Bundle::size() const {
  boost::unique_lock<boost::mutex>(_mut);
  if (_merged_into) {
    return _targets.size() + _merged_into->size();
  }
  return _targets.size();
}

void Bundle::incr_counts(size_t incr_amt) {
  boost::unique_lock<boost::mutex>(_mut);
  if (_merged_into) {
    _merged_into->incr_counts(incr_amt);
  } else {
    _counts += incr_amt;
  }
}

void Bundle::incr_mass(double incr_amt) {
  boost::unique_lock<boost::mutex>(_mut);
  if (_merged_into) {
    _merged_into->incr_mass(incr_amt);
  } else {
    _mass = log_add(_mass, incr_amt);
  }
}

void Bundle::reset_mass() {
  boost::unique_lock<boost::mutex>(_mut);
  _mass = LOG_0;
}

size_t Bundle::counts() const {
  boost::unique_lock<boost::mutex>(_mut);
  if (_merged_into) {
    return _merged_into->counts();
  }
  return _counts;
}

double Bundle::mass() const {
  boost::unique_lock<boost::mutex>(_mut);
  if (_merged_into) {
    return _merged_into->mass();
  }
  return _mass;
}


BundleTable::BundleTable() : _threadsafe_mode(false) {}

BundleTable::~BundleTable() {
  foreach(Bundle* bundle, _bundles) {
    delete bundle;
  }
}

Bundle* BundleTable::get_rep(Bundle* b) {
  if (b->_merged_into) {
    return get_rep(b->_merged_into);
  }
  return b;
}

Bundle* BundleTable::create_bundle(Target* targ) {
  Bundle* b = new Bundle(targ);
  _bundles.insert(b);
  return b;
}

Bundle* BundleTable::merge(Bundle* b1, Bundle* b2) {
  // Lock so that only one merge can happen at a time...
  boost::unique_lock<boost::mutex>(_mut);
  
  b1 = get_rep(b1);
  b2 = get_rep(b2);
  if (b1==b2) {
    return b1;
  }
  
  if (b1->size() < b2->size()) {
    swap(b1, b2);
  }
  
  if (_threadsafe_mode) {
    // Lock b1 and b2
    boost::unique_lock<boost::mutex>(b1->_mut);
    boost::unique_lock<boost::mutex>(b2->_mut);
    b1->_counts += b2->_counts;
    b1->_mass = log_add(b1->_mass, b2->_mass);
    b2->_counts = 0;
    b2->_mass = LOG_0;
    b2->_merged_into = b1;
  } else {
    foreach(Target* targ, b2->_targets) {
      targ->bundle(b1);
      b1->_targets.push_back(targ);
    }

    b1->incr_counts(b2->counts());
    b1->incr_mass(b2->mass());
    _bundles.erase(b2);
    delete b2;
  }
  
  return b1;
}

void BundleTable::collapse() {
  // Lock
  boost::unique_lock<boost::mutex>(_mut);
  
  BundleSet to_delete;
  foreach(Bundle* b, _bundles) {
    Bundle* rep = get_rep(b);
    if (rep != b) {
      foreach(Target* targ, b->_targets) {
        targ->bundle(rep);
        rep->_targets.push_back(targ);
      }
      to_delete.insert(b);
    }
  }
  
  foreach(Bundle* b, to_delete) {
    _bundles.erase(b);
    delete b;
  }
  
}

/**
 *  rhotree.h
 *  express
 *
 *  Created by Adam Roberts on 9/17/12.
 *
 **/

#ifndef __express__rhotree__
#define __express__rhotree__

#include <vector>

#include "boost/scoped_ptr.hpp"
#include "frequencymatrix.h"

class Fragment;

struct SapData {
  std::vector<size_t> leaf_ids;
  std::vector<double> rhos;
  std::vector<double> const_likelihoods;
  std::vector<double> accum_assignments;
  SapData(size_t size)
      : leaf_ids(size),
        rhos(size),
        const_likelihoods(size),
        accum_assignments(size + 1, LOG_0) {
  }
};

class Sap {
  SapData* _params;
  size_t _l;
  size_t _r;
  
 public:
  Sap (SapData* params, size_t l, size_t r) : _params(params), _l(l), _r(r) {}
  size_t size() const { return _r + 1 - _l; }
  size_t leaf_id(size_t i) const {
    return _params->leaf_ids[i];
  }
  double& rho(size_t i) {
    assert(i < size());
    return _params->rhos[i];
  }
  double const_likelihood(size_t i) const {
    assert(i < size());
    return _params->const_likelihoods[i];
  }
  double fraction() const { return log_sub(_params->accum_assignments[_r + 1],
                                           _params->accum_assignments[_l]); }
  Sap branch(size_t split);
};

class RhoTree {
  size_t _n;
  double _mass;
 protected:
  double _alpha;
  double _ff_param;
  std::vector<RhoTree*> _children;
  FrequencyMatrix<double> _child_rhos;
 public:
  RhoTree(double alpha, double ff_param)
      : _alpha(alpha), _n(0), _mass(0), _ff_param(ff_param) {}
  virtual ~RhoTree() {
    foreach(RhoTree* child, _children) {
      delete child;
    }
  }
  void fix() {
    _child_rhos = FrequencyMatrix<double>(1, _children.size(), _alpha);
  }
  double next_mass() {
    _n++;
    _mass += _ff_param*log((double)_n) - log(pow(_n+1,_ff_param) - 1);
    return _mass;
  }
  double similarity_scalar(const Sap& sap) { return 0; } // FIXME
  virtual void add_child(RhoTree* child) { _children.push_back(child); }
  virtual bool is_leaf() const = 0;
  virtual void get_rhos(Sap sap, double rho) const = 0;
  virtual void update_rhos(Sap sap) = 0;
};

class RangeRhoTree : public RhoTree {
  size_t _left;
  size_t _right;
 public:
  RangeRhoTree(size_t left, size_t right, double alpha, double ff_param);
  bool is_leaf() const { return _left == _right; }
  size_t left() const { return _left; }
  size_t right() const { return _right; }
  void add_child(RangeRhoTree* child);
  void get_rhos(Sap frag, double rho) const;
  void update_rhos(Sap frag);
};

typedef size_t LeafID;
typedef size_t TreeID;

class RhoForest : public RhoTree {
  // target-to-leaf (id)
  std::vector<LeafID> _target_to_leaf_map;
  // leaf-to-tree (id)
  std::vector<TreeID> _leaf_to_tree_map;
  void load_from_file(std::string infile);
 public:
  RhoForest(std::string infile, double alpha, double ff_param);
  void process_fragment(Fragment* frag);
  void get_rhos(Sap sap, double rho) { assert(false); }
  void update_rhos(Sap sap) { assert(false); }
};

#endif /* defined(__express__rhotree__) */

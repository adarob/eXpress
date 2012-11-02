/**
 *  rhotree.h
 *  express
 *
 *  Created by Adam Roberts on 9/17/12.
 *
 **/

#ifndef __express__rhotree__
#define __express__rhotree__

#include <stack>
#include <vector>

#include "boost/scoped_ptr.hpp"
#include "frequencymatrix.h"

class Fragment;

typedef size_t LeafID;
typedef size_t TreeID;

/**
 * Struct to hold values that are passed through the tree by the Sap class.
 */
struct SapData {
  TreeID tree_root;
  std::vector<LeafID> leaf_ids;
  std::vector<double> rhos;
  std::vector<double> const_likelihoods;
  std::vector<double> accum_const_likelihoods;
  std::vector<double> accum_assignments;
  SapData(size_t size)
      : leaf_ids(size),
        rhos(size),
        const_likelihoods(size),
        accum_const_likelihoods(size + 1, LOG_0),
        accum_assignments(size + 1, LOG_0) {
  }
};

/**
 * Class used to move data through the tree. Passed down to compute and update
 * rhos and likelihoods.
 */
class Sap {
  SapData* _params;
  size_t _l;
  size_t _r;

 public:
  Sap (SapData* params) : _params(params),
                          _l(0),
                          _r(_params->leaf_ids.size()-1){}
  Sap (SapData* params, size_t l, size_t r) : _params(params), _l(l), _r(r) {}
  TreeID tree_root() const { return _params->tree_root; }
  size_t size() const { return _r + 1 - _l; }
  size_t leaf_id(size_t i) const {
    assert(i < size());
    return _params->leaf_ids[_l + i];
  }
  double& rho(size_t i) {
    assert(i < size());
    return _params->rhos[_l + i];
  }
  double const_likelihood(size_t i) const {
    assert(i < size());
    return _params->const_likelihoods[_l + i];
  }
  double total_const_likelihood() const {
    return log_sub(_params->accum_const_likelihoods[_r + 1],
                   _params->accum_const_likelihoods[_l]);
  }
  double fraction() const {
    return log_sub(_params->accum_assignments[_r + 1],
                   _params->accum_assignments[_l]);
  }
  Sap branch(size_t split);
};

class RhoTree {
  size_t _n;
  double _mass;

protected:
  double _ff_param;
  std::vector<RhoTree*> _children;
  FrequencyMatrix<double> _child_rhos;

  virtual void add_child(RhoTree* child) { _children.push_back(child); }
  virtual void set_rhos(Sap sap) = 0;
  virtual void get_rhos(Sap sap, double rho) const = 0;
  virtual void update_rhos(Sap sap) = 0;
  double next_mass() {
    double ret = _mass;
    _n++;
    _mass += _ff_param*log((double)_n) - log(pow(_n+1,_ff_param) - 1);
    return ret;
  }
  double similarity_scalar(const Sap& sap);

public:
  friend class RhoLeafIterator;
  RhoTree(double ff_param)
      : _n(0), _mass(0), _ff_param(ff_param) {}
  virtual ~RhoTree() {
    foreach(RhoTree* child, _children) {
      delete child;
    }
  }
  size_t num_children() const { return _children.size(); }
  RhoTree* child(TreeID t) const { return _children[t]; };
  virtual size_t num_leaves() const = 0;
  virtual bool is_leaf() const = 0;
  virtual LeafID id() const = 0;
};

struct RhoLeafIteratorData {
  RhoTree* tree;
  size_t curr_child;
  double rho;
  RhoLeafIteratorData(RhoTree* t) : tree(t), curr_child(0), rho(0) {}
};

class RhoLeafIterator : public std::iterator<std::forward_iterator_tag,
                                             RhoTree> {
  std::stack<RhoLeafIteratorData> _stack;

  void prune_used_branches() {
    RhoLeafIteratorData* parent = NULL;
    while(parent == NULL) {
      if (_stack.empty()) {
        return;
      }
      parent = &_stack.top();
      parent->curr_child++;
      if (parent->curr_child >= parent->tree->num_children()) {
        _stack.pop();
        parent = NULL;
      }
    }
  }
  void travel_to_leaf() {
    RhoLeafIteratorData* parent = &_stack.top();
    while (!parent->tree->is_leaf()) {
      RhoLeafIteratorData child(parent->tree->_children[parent->curr_child]);
      child.rho = parent->rho + parent->tree->_child_rhos(parent->curr_child);
      _stack.push(child);
      parent = &_stack.top();
    }
  }

 public:
  RhoLeafIterator(RhoTree* tree) {
    RhoLeafIteratorData d(tree);
    _stack.push(d);
    travel_to_leaf();
  }
  RhoLeafIterator& operator++() {
    prune_used_branches();
    if (!_stack.empty()) {
      travel_to_leaf();
    }
    return *this;
  }
  RhoTree* operator*() {
    if (_stack.empty()) {
      return NULL;
    }
    return _stack.top().tree;
  }
  double rho() {
    if (_stack.empty()) {
      return LOG_0;
    }
    return _stack.top().rho;
  }
};

class RangeRhoTree : public RhoTree {
 protected:
  friend class RangeRhoForest;
  size_t _left;
  size_t _right;
  void add_child(RangeRhoTree* child);
  void set_rhos(Sap sap);
  virtual void get_rhos(Sap sap, double rho) const;
  virtual void update_rhos(Sap sap);
 public:
  RangeRhoTree(size_t left, size_t right, double ff_param);
  bool is_leaf() const { return _left == _right; }
  size_t left() const { return _left; }
  size_t right() const { return _right; }
  size_t num_leaves() const { return _right - _left + 1; }
  LeafID id() const {
    if (is_leaf()) {
      return _left;
    }
    return -1;
  }
};

class RangeRhoForest : public RangeRhoTree {
  // target-to-leaf (id)
  std::vector<LeafID> _target_to_leaf_map;
  // leaf-to-tree (id)
  std::vector<TreeID> _leaf_to_tree_map;
  std::vector<size_t> _tree_counts;
  void load_from_file(std::string infile);
protected:
  void get_rhos(Sap sap, double rho) const;
  void update_rhos(Sap sap);
public:
  RangeRhoForest(std::string infile, double ff_param);
  void set_alphas(const std::vector<double>& target_alphas);
  void process_fragment(const Fragment& frag);
  size_t tree_counts(TreeID t) const { return _tree_counts[t]; }
  size_t num_leaves() const {
    return _leaf_to_tree_map.size();
  }
  // TODO: We stop using this map after remapping the targets...need to fix this.
  const std::vector<LeafID>* target_to_leaf_map() const {
    return &_target_to_leaf_map;
  }
};

#endif /* defined(__express__rhotree__) */

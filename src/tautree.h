/**
 *  tautree.h
 *  express
 *
 *  Created by Adam Roberts on 9/17/12.
 *
 **/

#ifndef __express__tautree__
#define __express__tautree__

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
  std::vector<double> taus;
  std::vector<double> const_likelihoods;
  std::vector<double> accum_const_likelihoods;
  std::vector<double> accum_assignments;
  SapData() {}
  SapData(size_t size)
      : leaf_ids(size),
        taus(size),
        const_likelihoods(size),
        accum_const_likelihoods(size + 1, LOG_0),
        accum_assignments(size + 1, LOG_0) {
  }
};

/**
 * Class used to move data through the tree. Passed down to compute and update
 * taus and likelihoods.
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
  double& tau(size_t i) {
    assert(i < size());
    return _params->taus[_l + i];
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

class TauTree {
  // logged
  double _n;
  double _mass;

protected:
  double _ff_param;
  std::vector<TauTree*> _children;
  FrequencyMatrix<double> _child_taus;

  virtual void add_child(TauTree* child) { _children.push_back(child); }
  virtual void get_taus(Sap sap, double tau) const = 0;
  virtual void update_taus(Sap sap) = 0;
  virtual void increment_taus(Sap sap, bool decrement=false) = 0;
  double next_mass(double incr) {
    double new_n = log_add(_n, incr);
    _mass += _ff_param*_n - log_sub(_ff_param*new_n, LOG_1);
    _n = new_n;
    return _mass;
  }
  double similarity_scalar(const Sap& sap);

public:
  friend class TauLeafIterator;
  TauTree(double ff_param)
      : _n(LOG_1), _mass(LOG_1), _ff_param(ff_param) {}
  virtual ~TauTree() {
    foreach(TauTree* child, _children) {
      delete child;
    }
  }
  size_t num_children() const { return _children.size(); }
  TauTree* child(TreeID t) const { return _children[t]; };
  virtual size_t num_leaves() const = 0;
  virtual bool is_leaf() const = 0;
  virtual LeafID id() const = 0;
  virtual void renormalize_taus();
};

struct TauLeafIteratorData {
  TauTree* tree;
  size_t curr_child;
  double tau;
  TauLeafIteratorData(TauTree* t) : tree(t), curr_child(0), tau(0) {}
};

class TauLeafIterator : public std::iterator<std::forward_iterator_tag,
                                             TauTree> {
  std::stack<TauLeafIteratorData> _stack;

  void prune_used_branches() {
    TauLeafIteratorData* parent = NULL;
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
    TauLeafIteratorData* parent = &_stack.top();
    while (!parent->tree->is_leaf()) {
      TauLeafIteratorData child(parent->tree->_children[parent->curr_child]);
      child.tau = parent->tau + parent->tree->_child_taus(parent->curr_child);
      _stack.push(child);
      parent = &_stack.top();
    }
  }

 public:
  TauLeafIterator(TauTree* tree) {
    TauLeafIteratorData d(tree);
    _stack.push(d);
    travel_to_leaf();
  }
  TauLeafIterator& operator++() {
    prune_used_branches();
    if (!_stack.empty()) {
      travel_to_leaf();
    }
    return *this;
  }
  TauTree* operator*() {
    if (_stack.empty()) {
      return NULL;
    }
    return _stack.top().tree;
  }
  double tau() {
    if (_stack.empty()) {
      return LOG_0;
    }
    return _stack.top().tau;
  }
};

class RangeTauTree : public TauTree {
 protected:
  friend class RangeTauForest;
  size_t _left;
  size_t _right;
  void add_child(RangeTauTree* child);
  void initialize_taus();
  virtual void get_taus(Sap sap, double tau) const;
  virtual void update_taus(Sap sap);
  virtual void increment_taus(Sap sap, bool decrement=false);
 public:
  RangeTauTree(size_t left, size_t right, double ff_param);
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

class RangeTauForest : public RangeTauTree {
  // target-to-leaf (id)
  std::vector<LeafID> _target_to_leaf_map;
  // leaf-to-tree (id)
  std::vector<TreeID> _leaf_to_tree_map;
  std::vector<size_t> _tree_counts;
  SapData _alpha_params;
  void load_from_file(std::string infile);
protected:
  void get_taus(Sap sap, double tau) const;
  void update_taus(Sap sap);
public:
  RangeTauForest(std::string infile, double ff_param);
  void set_alphas(const std::vector<double>& target_alphas);
  void add_alphas() { increment_taus(Sap(&_alpha_params)); };
  void remove_alphas() {
    increment_taus(Sap(&_alpha_params), true);
    renormalize_taus();
  }
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

#endif /* defined(__express__tautree__) */

/**
 *  tautree.cpp
 *  express
 *
 *  Created by Adam Roberts on 9/17/12.
 *
 **/

#include "tautree.h"

#include "main.h"
#include "fragments.h"
#include "targets.h"

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

using namespace std;

Sap Sap::branch(size_t split) {
  size_t old_left = _l;
  vector<size_t>::iterator it = upper_bound(_params->leaf_ids.begin() + _l,
                                            _params->leaf_ids.begin() + _r + 1,
                                            split);
  _l = it - _params->leaf_ids.begin();
  return Sap(_params, old_left, _l-1);
}

double TauTree::similarity_scalar(const Sap& sap) {
  if (sap.size() < _children.size()) {
    return LOG_1;
  }
  double c = 0;
  double tot = sap.total_const_likelihood();
  if (tot == LOG_0) {
    return LOG_EPSILON;
  }
  for (size_t i = 0; i < sap.size(); ++i) {
    double p = sap.const_likelihood(i) - tot;
    c += -sexp(p)*p;
    assert(!isnan(c));
  }
  c = 1 - (c / log(sap.size()));
  if (c < 0 || approx_eq(c,0)) {
    return LOG_EPSILON;
  }
  if (c > 1) {
    return LOG_1;
  }
  assert(c==0 || !isnan(log(c)));
  return log(c);
}

void TauTree::renormalize_taus() {
  _child_taus.recompute_sums();
  foreach(TauTree* child, _children) {
    child->renormalize_taus();
  }
}

RangeTauForest::RangeTauForest(string infile_path, double ff_param)
    : RangeTauTree(0, 0, ff_param) {
  load_from_file(infile_path);
}

void RangeTauForest::load_from_file(string infile_path) {
  ifstream infile (infile_path.c_str());
  string line;

  if (!infile.is_open()) {
    cerr << "ERROR: Could not open hierarchy specification file '"
         << infile_path << "'.";
    exit(1);
  }

  getline(infile, line, '\n');
  size_t split = line.find(',');
  size_t num_leaves = atoi(line.substr(0, split).c_str());
  size_t num_nodes = atoi(line.substr(split+1).c_str());
  _target_to_leaf_map = vector<LeafID>(num_leaves, -1);
  _leaf_to_tree_map = vector<TreeID>(num_leaves, -1);
  LeafID next_leaf_id = 0;
  vector<RangeTauTree*> nodes(num_nodes, NULL);

  cout << "Loading hierarchy from '" << infile_path << "' with " << num_nodes
       << " nodes and " << num_leaves << " leaves...\n";
  boost::char_separator<char> sep_edges(";");

  while (infile.good()) {
    getline(infile, line, '\n');
    boost::algorithm::trim(line);
    if (line.empty()) {
      continue;
    }
    boost::tokenizer<boost::char_separator<char> > tokens(line, sep_edges);
    TreeID parent, child;
    foreach (const string& edge, tokens) {
      size_t split = edge.find(',');
      parent = atoi(edge.substr(0, split).c_str());
      child = atoi(edge.substr(split+1).c_str());
      if (child < num_leaves) {
        assert(_target_to_leaf_map[child] == -1);
        nodes[next_leaf_id] = new RangeTauTree(next_leaf_id, next_leaf_id,
                                               _ff_param);
        _target_to_leaf_map[child] = next_leaf_id;
        _leaf_to_tree_map[next_leaf_id] = _children.size();
        child = next_leaf_id;
        next_leaf_id++;
      }
      assert(nodes[child]);
      if (!nodes[parent]) {
        nodes[parent] = new RangeTauTree(nodes[child]->left(),
                                         nodes[child]->right(),
                                         _ff_param);
      }
      else if (nodes[parent]->right() + 1 != nodes[child]->left()) {
        cerr << "ERROR: '" << infile_path << "' does not specify a proper "
             << "hierarchy.";
        exit(1);
      }
      nodes[parent]->add_child(nodes[child]);
    }
    add_child(nodes[parent]);
  }
  
  for(TargID t = 0; t < num_leaves; ++t) {
    if (_target_to_leaf_map[t] == -1) {
      nodes[next_leaf_id] = new RangeTauTree(next_leaf_id, next_leaf_id,
                                             _ff_param);
      _target_to_leaf_map[t] = next_leaf_id;
      _leaf_to_tree_map[next_leaf_id] = _children.size();
      add_child(nodes[next_leaf_id]);
      next_leaf_id++;
    }
  }
  
  assert(next_leaf_id == num_leaves);
  _right = num_leaves-1;
  _tree_counts = vector<size_t>(_children.size(), 0);
  initialize_taus();
}

void RangeTauForest::set_alphas(const vector<double>& target_alphas) {
  _alpha_params = SapData(target_alphas.size());
  assert(target_alphas.size() == num_leaves());
  
  for (TargID i = 0; i < target_alphas.size(); i++) {
    //    LeafID leaf = _target_to_leaf_map[i];
    LeafID leaf = i;
    _alpha_params.leaf_ids[leaf] = leaf;
    _alpha_params.taus[leaf] = target_alphas[i];
    assert(_alpha_params.taus[leaf] != LOG_0);
  }
  for (LeafID i = 0; i < num_leaves(); ++i) {
    assert(_alpha_params.taus[i] != LOG_0);
    _alpha_params.accum_assignments[i+1] = log_add(_alpha_params.accum_assignments[i],
                                                   _alpha_params.taus[i]);
  }
  add_alphas();
}

void RangeTauForest::get_taus(Sap sap, double tau) const {
  TreeID tree = sap.tree_root();
  static_cast<RangeTauTree*>(_children[tree])->get_taus(sap, _child_taus(tree));
}

void RangeTauForest::update_taus(Sap sap) {
  TreeID tree = sap.tree_root();

  // Don't compute a scalar here to save time. These should always be valuable.
  double mass = next_mass();
  // Update the taus throughout the forest using the computed likelihoods.
  assert(approx_eq(sap.fraction(), LOG_1));
  _child_taus.increment(tree, mass);
  static_cast<RangeTauTree*>(_children[tree])->update_taus(sap);
}

void RangeTauForest::process_fragment(const Fragment& frag) {
  // Assume the likelihoods are pre-computed in probability field.

  // Everything is much easier and faster when the mapping is unique.
  if (frag.num_hits() == 1) {
    SapData params(1);
    FragHit& hit = *frag[0];
    hit.probability = LOG_1;
    params.leaf_ids[0] = hit.targ_id;
    params.tree_root = _leaf_to_tree_map[params.leaf_ids[0]];
    params.accum_assignments[1] = LOG_1;
    double mass = next_mass();
    _child_taus.increment(params.tree_root, mass);
    static_cast<RangeTauTree*>(_children[params.tree_root])->update_taus(Sap(&params, 0, 0));
    _tree_counts[params.tree_root]++;
    return;
  }

  SapData params(frag.num_hits());

  for (size_t i = 0; i < frag.num_hits(); ++i) {
    FragHit& hit = *frag[i];
    params.leaf_ids[i] = hit.targ_id;
    params.const_likelihoods[i] = hit.probability;
    params.accum_const_likelihoods[i+1] = log_add(params.accum_const_likelihoods[i],
                                                 hit.probability);
    if (i == 0) {
      params.tree_root = _leaf_to_tree_map[params.leaf_ids[i]];
    } else {
      assert(_leaf_to_tree_map[params.leaf_ids[i]] == params.tree_root);
      if (_leaf_to_tree_map[params.leaf_ids[i]] != params.tree_root) {
        cerr << "ERROR: Invalid hierarchy forest. Alignment of '"
             << frag.name() << "' accesses multiple trees.";
        exit(1);
      }
    }
  }

  Sap sap(&params);
  get_taus(sap, LOG_1);

  // Compute the individual and total likelihoods.
  double total_likelihood = LOG_0;
  for (size_t i = 0; i < frag.num_hits(); ++i) {
    total_likelihood = log_add(total_likelihood,
                               params.const_likelihoods[i] + params.taus[i]);
  }

  // Accumulate the fractional assignments in the Sap for easy lookup later on.
  for (size_t i = 0; i < frag.num_hits(); ++i) {
    double frac = params.const_likelihoods[i] + params.taus[i] -
                  total_likelihood;
    frag[i]->probability = frac;
    params.accum_assignments[i+1] = log_add(params.accum_assignments[i], frac);
  }
  
  update_taus(sap);
  _tree_counts[params.tree_root]++;
}

RangeTauTree::RangeTauTree(size_t left, size_t right, double ff_param)
    : TauTree(ff_param), _left(left), _right(right) {
}

void RangeTauTree::add_child(RangeTauTree* child) {
  if (_children.empty()) {
    _left = child->left();
    _right = child->right();
  } else {
    assert(child->left() == _right + 1);
    _right = child->right();
  }
  TauTree::add_child(child);
}

void RangeTauTree::initialize_taus() {
  if (is_leaf()) {
    return;
  }
  _child_taus = FrequencyMatrix<double>(1, _children.size(), LOG_0);
  foreach (TauTree* child, _children) {
    static_cast<RangeTauTree*>(child)->initialize_taus();
  }
}

void RangeTauTree::get_taus(Sap sap, double tau) const {
  assert(!isnan(tau));
  if (is_leaf()) {
    for (size_t i = 0; i < sap.size(); ++i) {
      assert(sap.leaf_id(i) == id());
      sap.tau(i) = tau;
    }
    return;
  }

  for (size_t i = 0; i < _children.size(); ++i) {
    RangeTauTree& child = *static_cast<RangeTauTree*>(_children[i]);
    Sap branch_sap = sap.branch(child.right());
    if (branch_sap.size()) {
      child.get_taus(branch_sap, tau + _child_taus(i));
    }
    if (sap.size() == 0) {
      break;
    }
  }
}

void RangeTauTree::update_taus(Sap sap) {
  if (is_leaf()) {
    return;
  }
  double sim_scalar = similarity_scalar(sap);
  if (islzero(sim_scalar)) {
    return;
  }
  const double mass = sim_scalar + next_mass();
  
  for (size_t i = 0; i < _children.size(); ++i) {
    RangeTauTree& child = *static_cast<RangeTauTree*>(_children[i]);
    Sap branch_sap = sap.branch(child.right());
    if (branch_sap.size()) {
      child.update_taus(branch_sap);
      _child_taus.increment(i, mass + branch_sap.fraction());
    }
    if (sap.size() == 0) {
      break;
    }
  }
}

void RangeTauTree::increment_taus(Sap sap, bool decrement) {
  if (is_leaf()) {
    return;
  }
  for (size_t i = 0; i < _children.size(); ++i) {
    RangeTauTree& child = *static_cast<RangeTauTree*>(_children[i]);
    Sap branch_sap = sap.branch(child.right());
    if (branch_sap.size()) {
      child.increment_taus(branch_sap, decrement);
      if (decrement) {
        _child_taus.decrement(i, branch_sap.fraction());
      } else {
        _child_taus.increment(i, branch_sap.fraction());
      }
    }
    if (sap.size() == 0) {
      break;
    }
  }
}

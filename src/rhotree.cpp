/**
 *  rhotree.cpp
 *  express
 *
 *  Created by Adam Roberts on 9/17/12.
 *
 **/

#include "rhotree.h"

#include "main.h"
#include "fragments.h"
#include "targets.h"

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

void RangeRhoForest::load_from_file(string infile_path) {
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
  LeafID next_leaf_id = 0;
  vector<RangeRhoTree*> nodes(num_nodes, NULL);
  
  cout << "Loading hierarchy from '" << infile_path << "' with " << num_nodes
       << " nodes and " << num_leaves << "leaves...\n";
  boost::char_separator<char> sep_edges(";");
  
  while (infile.good()) {
    getline(infile, line, '\n');
    boost::tokenizer<boost::char_separator<char> > tokens(line, sep_edges);
    TreeID parent, child;
    foreach (const string& edge, tokens) {
      size_t split = edge.find(',');
      parent = atoi(edge.substr(0, split).c_str());
      child = atoi(edge.substr(split+1).c_str());
      if (child < num_leaves) {
        assert(!_target_to_leaf_map[child]);
        nodes[next_leaf_id] = new RangeRhoTree(child, child, _ff_param);
        _target_to_leaf_map[child] = next_leaf_id;
        _leaf_to_tree_map[next_leaf_id] = _children.size();
        child = next_leaf_id;
        next_leaf_id++;
      }
      assert(nodes[child]);
      if (!nodes[parent]) {
        nodes[parent] = new RangeRhoTree(nodes[child]->left(),
                                         nodes[child]->right(),
                                         _ff_param);
      }
      if (nodes[parent]->right() != nodes[child]->left() + 1) {
        cerr << "ERROR: '" << infile_path << "' does not specify a proper "
             << "hierarchy.";
        exit(1);
      }
      nodes[parent]->add_child(nodes[child]);
    }
    add_child(nodes[parent]);
  }
}

SapData RangeRhoForest::initialize_sapdata(Fragment* frag) {
  SapData params(frag->num_hits());
  
  for (size_t i = 0; i < frag->num_hits()) {
    FragHit& hit = *(frag->hits()[i]);
    params.leaf_ids[i] = _target_to_leaf_map[hit.targ_id]
    if (i == 0) {
      params.tree = _leaf_to_tree_map[_params.leaf_ids[i]];
    } else {
      assert(_leaf_to_tree_map[_params.leaf_ids[i]] == params.tree);
    }
  }
}

void RangeRhoForest::set_alphas(vector<double> target_alphas) {
  SapData params(target_alphas.size());
  for (size_t i = 0; i < target_alphas.size(); i++) {
    LeafID leaf = _target_to_leaf_map[i];
    params.leaf_ids[leaf] = leaf;
    params.rhos[leaf] = target_alphas[i];
  }
  set_rhos(Sap(&params));
}

void RangeRhoForest::get_rhos(Sap sap) {
  TreeID tree = sap.tree_root();
  _children[tree]->get_rhos(sap, _child_rhos(tree));
  return params.rhos;
}

void RangeRhoForest::update_rhos(Sap sap) {
  TreeID tree = sap.tree_root();
  
  // Don't compute a scalar here to save time. These should always be valuable.
  double mass = next_mass();
  // Update the rhos throughout the forest using the computed likelihoods.
  _child_rhos.increment(tree_id, mass + tree_saps[j].fraction());
  _children[tree]->update_rhos(tree_saps[j]);
}

void RangeRhoForest::process_fragment(Fragment* frag) {
  // TODO: Assume the likelihoods are pre-computed, perhaps pass in through the probabilities field?
  
  // Everything is much easier and faster when the mapping is unique.
  if (frag->num_hits() == 1) {
    SapData params(1);
    FragHit& hit = *frag->hits()[0];
    params.leaf_ids[0] = _target_to_leaf_map[hit.targ_id]
    params.tree = _leaf_to_tree_map[params.leaf_ids[0]]];
    params.accum_assignments[1] = 0;
    double mass = next_mass();
    _child_rhos.increment(params.tree, mass);
    _children[tree_id]->update_rhos(Sap(&params, 0, 0));
    return;
  }

  SapData params = initialize_sapdata(frag);
  // TODO: sort leaf_ids and create a mapping to fraghits
  Sap sap(&params);
  get_rhos(sap);
  
  // Compute the individual and total likelihoods.
  // TODO: change this to use the mapping and the probability field
  double total_likelihood = LOG_0;
  for (size_t i = 0; i < frag->num_hits(); ++i) {
    FragHit& hit = *frag->hits()[i];
    params.const_likelihoods[i] = hit.targ->log_likelihood(hit, true);
    total_likelihood = log_add(total_likelihood,
                               params.const_likelihoods[i] + params.rhos[i]);
  }
  
  // Accumulate the fractional assignments in the Sap for easy lookup later on.
  for (size_t i = 0; i < frag->num_hits(); ++i) {
    double frac = params.const_likelihoods[i] + params.rhos[i] -
                  total_likelihood;
    frag[i]->probability = frac;
    params.accum_assignments[i+1] = log_add(frac, params.accum_assignments[i]);
  }
  
  update_rhos(sap);
}

RangeRhoTree::RangeRhoTree(size_t left, size_t right, double ff_param)
    : RhoTree(ff_param), _left(left), _right(right) {
}

void RangeRhoTree::add_child(RangeRhoTree* child) {
  if (_children.empty()) {
    _left = child->left();
    _right = child->right();
  } else {
    assert(child->left() == _right + 1);
    _right = child->right();
  }
  RhoTree::add_child(child);
}

RangeRhoTree::set_rhos(Sap sap) {
  _child_rhos = FrequencyMatrix<double>(1, _children.size(), LOG_0);
  for (size_t i = 0; i < _children.size(); ++i) {
    RangeRhoTree& child = *static_cast<RangeRhoTree*>(_children[i]);
    Sap branch_sap = sap.branch(child.right());
    if (branch_sap.size()) {
      child.set_rhos(branch_sap);
      _child_rhos.increment(i, branch_sap.rhos(i));
    }
    if (sap.size() == 0) {
      break;
    }
  }
}

void RangeRhoTree::get_rhos(Sap sap, double rho) const {
  if (is_leaf()) {
    for (size_t i = 0; i < sap.size(); ++i) {
      sap.rho(i) = rho;
    }
    return;
  }

  for (size_t i = 0; i < _children.size(); ++i) {
    RangeRhoTree& child = *static_cast<RangeRhoTree*>(_children[i]);
    Sap branch_sap = sap.branch(child.right());
    if (branch_sap.size()) {
      child.get_rhos(branch_sap, rho + _child_rhos(i));
    }
    if (sap.size() == 0) {
      break;
    }
  }

}

void RangeRhoTree::update_rhos(Sap sap) {
  const double mass = next_mass() + similarity_scalar(sap);
  
  for (size_t i = 0; i < _children.size(); ++i) {
    RangeRhoTree& child = *static_cast<RangeRhoTree*>(_children[i]);
    Sap branch_sap = sap.branch(child.right());
    if (branch_sap.size()) {
      child.update_rhos(branch_sap);
      _child_rhos.increment(i, mass + branch_sap.fraction());
    }
    if (sap.size() == 0) {
      break;
    }
  }
}
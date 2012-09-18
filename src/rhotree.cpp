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

using namespace std;

Sap Sap::branch(size_t split) {
  size_t old_left = _l;
  vector<size_t>::iterator it = upper_bound(_params->leaf_ids.begin() + _l,
                                            _params->leaf_ids.begin() + _r + 1,
                                            split);
  _l = it - _params->leaf_ids.begin();
  return Sap(_params, old_left, _l-1);
}


void RhoForest::process_fragment(Fragment* frag) {
  SapData params(frag->num_hits());
  
  // Everything is much easier and faster when the mapping is unique.
  if (frag->num_hits() == 1) {
    FragHit& hit = *frag->hits()[0];
    hit.targ_id = _leaf_map[hit.targ_id];
    params.leaf_ids[0] = hit.targ_id;
    TreeID tree_id = _tree_map[hit.targ_id];
    params.accum_assignments[1] = 0;
    double mass = next_mass();
    _child_rhos.increment(tree_id, mass);
    _children[tree_id]->update_rhos(Sap(&params, 0, 0));
    return;
  }

  // These are the relevant roots and associated saps.
  vector<TreeID> tree_roots;
  vector<Sap> tree_saps;
  
  const TreeID NONE = _children.size() + 1;
  TreeID curr_tree = NONE;
  
  // Look up leaf ids and find the relevant trees.
  size_t next_left = 0;
  for (size_t i = 0; i < frag->num_hits(); ++i) {
    FragHit& hit = *frag->hits()[i];
    // This is not very elegant since it replaces the target id with the leaf
    // id, but it allows for easy sorting.
    hit.targ_id = _leaf_map[hit.targ_id];
    params.leaf_ids[i] = hit.targ_id;
    TreeID tree = _tree_map[hit.targ_id];
    if (tree != curr_tree) {
      if (curr_tree != NONE) {
        tree_roots.push_back(curr_tree);
        tree_saps.push_back(Sap(&params, next_left, i));
        next_left = i + 1;
      }
      curr_tree = tree;
    }
  }
  tree_roots.push_back(curr_tree);
  tree_saps.push_back(Sap(&params, next_left, frag->num_hits()-1));
  
  // For now we'll assume the full trees are bundles, although nothing except
  // this assertions relies on that assumption.
  assert(tree_roots.size() == 0);

  // Use breadth-first search to get the rhos
  for (size_t j = 0; j < tree_roots.size(); ++j) {
    size_t tree_id = tree_roots[j];
    _children[tree_id]->get_rhos(tree_saps[j], _child_rhos(tree_id));
  }
  
  // Sort the FragHits based on leaf_id.
  frag->sort_hits();
  
  // Compute the individual and total likelihoods.
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
    params.accum_assignments[i+1] = log_add(frac, params.accum_assignments[i]);
  }
    
  // Don't compute a scalar here to save time. These should always be valuable.
  double mass = next_mass();
  // Update the rhos throughout the forest using the computed likelihoods.
  for (size_t j = 0 ; j < tree_roots.size(); ++j) {
    size_t tree_id = tree_roots[j];
    _child_rhos.increment(tree_id, mass + tree_saps[j].fraction());
    _children[tree_id]->update_rhos(tree_saps[j]);
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
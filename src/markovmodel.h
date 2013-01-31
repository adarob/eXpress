/**
 *  markovmodel.h
 *  express
 *
 *  Created by Adam Roberts on 1/24/12.
 *  Copyright 2012 Adam Roberts. All rights reserved.
 **/

#ifndef express_markovmodel_h
#define express_markovmodel_h

#include <vector>
#include <string>
#include "frequencymatrix.h"

class Sequence;

/**
 * The MarkovModel class is used to store transition probabilities of a Markov
 * chain based on a nucleotide Sequence, which itself can be probabilistic. The
 * The probabilities can be updated based on a sequence window and can be
 * initialized by sliding a window over a sequence. Nucleotides are represented
 * by size_t values according to the NUCS array in main.h.
 *  @author  Adam Roberts
 *  @date    2012
 *  @copyright Artistic License 2.0
 **/
class MarkovModel {
  /**
   * A private int storing the order of the Markov model.
   */
  int _order;
  /**
   * A private int storing the size of the sequence window.
   */
  int _window_size;
  /**
   * A private int storing the number of nodes in the model.
   */
  int _num_pos;
  /**
   * A vector of FrequencyMatrix objects to store the transition probabilities.
   */
  std::vector<FrequencyMatrix<double> > _params;
  /**
   * A size_t used to clear higher-order bits when iterating through conditional
   * probabilities indices.
   */
  size_t _bitclear;

 public:
  /**
   * MarkovModel Constructor.
   * @param order the order of the model.
   * @param window_size the size of the sequence window to calculate
   *        probabilities for.
   * @param num_pos the number of nodes in the model.
   * @param alpha the initial pseudo-counts (non-logged).
   */
  MarkovModel(size_t order, size_t window_size, size_t num_pos, double alpha);
  /**
   * Accessor for the probability of transitioning from cond to curr at node
   * p.
   * @param p the node to get the transition probability from.
   * @param cond the index of the previous state.
   * @param curr the index of the state being transitioned to.
   * @return The probability of transitioning from cond to curr at node p.
   */
  double transition_prob(size_t p, size_t cond, size_t curr) const;
  /**
   * Computes the probability of the sequence beginning at left of size
   * _window_size using the parameters of the Markov model.
   * @param seq the sequence from which to extract the window.
   * @param left the leftmost point in the sequence window.
   * @return The probability of the sequence based on the model parameters.
   */
  double seq_prob(const Sequence& seq, int left) const;
  /**
   * Increments the parameters associated with the sequence beginning at left of
   * size _window_size by the (logged) mass.
   * @param seq the sequence from which to extract the window.
   * @param left the leftmost point in the sequence window.
   * @param mass the amount to increment the parameters by (logged).
   */
  void update(const Sequence& seq, int left, double mass);
  /**
   * Increments the specified transition by the given mass.
   * @param p the position in the chain to increment.
   * @param i the index of the previous state.
   * @param j the index of the transitioned state.
   * @param mass the amount to increment by (logged).
   */
  void update(size_t p, size_t i, size_t j, double mass);
  /**
   * A member function that computes and returns the parameter table indices
   * used to compute and update the likelihood for the given sequence at the
   * window started at the given position. Negative values indicate that the
   * position is not used in the likelihood computation.
   * @param seq the sequence to find the indices for.
   * @param left the position where the window should begin in the sequence.
   * @param indices a vector to fill with indices into the parameter tables that
   *        would be used in a likelihood calculation.
   * @return The starting position in the sequence to align with the first
   *         index.
   */
  size_t get_indices(const Sequence& seq, int left, std::vector<char>& indices);
  /**
   * A member function that computes and returns the parameter table indices
   * used to fast_learn (see below) the likelihood for the given sequence.
   * Negative values indicate that the position is not used in the fast_learn
   * computation.
   * @param seq the sequence to find the indices for.
   * @return A vector of indices into the parameter tables that would be used
   *         in a fast_learn.
   */
  std::vector<char> get_indices(const Sequence& seq);
  /**
   * Computes the marginal probability of transitioning to the given nucleotide
   * at position w in the model.
   * @param w the position (node) in the model.
   * @param nuc the nucleotide to calculate the marginal transition probability
   *        to.
   * @return The marginal probability of transitioning to nuc at position w.
   */
  double marginal_prob(size_t w, size_t nuc) const;
  /**
   * Slides a window along the given sequence, incrementing the highest order
   * transition paramaters by the given mass multiplied by the probability of
   * observing a fragment at that distance from the end.
   * @param seq the sequence to slide the window along.
   * @param mass the amount to increment the parameters by.
   * @param fl_cmf the fragment length CMF to determine the probability of
   *        observing fragment starts at different positions in the sequence.
   */
  void fast_learn(const Sequence& seq, double mass,
                  const std::vector<double>& fl_cmf);
  /**
   * After learning the highest order transitions with fast_learn, this method
   * fills in the lower-order transitions.
   */
  void calc_marginals();
};

#endif

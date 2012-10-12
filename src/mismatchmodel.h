#ifndef MISMATCH_H
#define MISMATCH_H

/**
 *  mismatchmodel.h
 *  express
 *
 *  Created by Adam Roberts on 4/23/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 **/

#include "frequencymatrix.h"

class FragHit;
class Target;

/**

 * The MismatchTable class is used to store and update mismatch and indel
 * (error) parameters using a first-order Markov model based on nucleotide and
 * position in a read. Also computes likelihoods of mismatches and indels in
 * given fragment mappings. When the target sequences are proabilistic, this
 * class is responsible for updating those parameters. All values are stored and
 * returned in log space.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class MismatchTable {
  /**
   * A vector of FrequencyMatrix objects to store the Markov model parameters
   * for each position in the first ("left") read.
   */
  std::vector<FrequencyMatrix<double> > _first_read_mm;
  /**
   * A vector of FrequencyMatrix objects to store the Markov model parameters
   * for each position in the second ("right") read.
   */
  std::vector<FrequencyMatrix<double> > _second_read_mm;
  /**
   * A FrequencyMatrix storing the observations of insertions of given lengths.
   */
  FrequencyMatrix<double> _insert_params;
  /**
   * A FrequencyMatrix storing the observations of deletions of given lengths.
   */
  FrequencyMatrix<double> _delete_params;
  /**
   * A size_t storing the maximum observed read length.
   */
  size_t _max_len;
  /**
   * A boolean specifying whether or not the table values are burned in. When
   * inactive, a log-likelihood is not computed (0 is returned) and
   * probabalistic target sequences are not updated.
   */
  bool _active;

 public:
  /**
   * MismatchTable constructor initializes the model parameters using the
   * specified (non-logged) pseudo-counts.
   * @param alpha a double containing the non-logged pseudo-counts for parameter
   *        initialization.
   */
  MismatchTable(double alpha);
  /**
   * Mutator to set the _active member variable to allow for log_likelihood
   * calculations. Used to skip calculations before burn-in completes.
   * @param active a boolean specifying whether to activate (true) or deactivate
   *        (false)
   */
  void activate(bool active = true) { _active = active; }
  /**
   * A member function that computes and returns the parameter table indices
   * used to compute and update the likelihood based on the given FragHit.
   * Negative values indicate that the position is not used in the likelihood
   * computation.
   * @param f the FragHit to find the table indices for.
   * @param left_indices the vector to store the left read's indices in.
   * @param right_indices the vector to store the right read's indices in.
   */
  void get_indices(const FragHit& f,
                   std::vector<char>& left_indices,
                   std::vector<char>& right_indices) const;
  /**
   * A member function that returns the log likelihood of mismatches and indels
   * in the mapping given the current error model parematers. Returns 0 if
   * _active is false.
   * @param f the fragment mapping to calculate the log likelihood for.
   * @return The log likelihood of the mapping based on mismatches and indels.
   */
  double log_likelihood(const FragHit& f) const;
  /**
   * A member function that updates the error model parameters based on a
   * mapping and its (logged) mass. Also updates the sequence parameters if
   * they are probabilistic and active_ is true.
   * @param f the fragment mapping.
   * @param p the logged posterior probablity of the alignment.
   * @param mass the logged mass of the fragment.
   */
  void update(const FragHit&, double p, double mass);
  /**
   * Freezes the parameters to allow for faster computation after burn out.
   * Cannot be undone.
   */
  void fix();
  /**
   * A member function that appends the final model parameters in tab-separated
   * format to the given file. The output has 1 row for each read position and
   * the parameters are in columns indexed as (ref, prev, obs) in base 4 with
   * A,C,G,T encoded as 0,1,2,3.
   * @param file stream to append to.
   */
  void append_output(std::ofstream& outfile) const;
};

#endif

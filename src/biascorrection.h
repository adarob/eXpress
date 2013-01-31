/**
 *  biascorrection.h
 *  express
 *
 *  Created by Adam Roberts on 4/5/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 *
 */

#ifndef BIASCORRECTION_H
#define BIASCORRECTION_H

#include <string>
#include <vector>

#include "lengthdistribution.h"
#include "frequencymatrix.h"
#include "markovmodel.h"

class BiasBoss;
class FragHit;
class Target;

/**
 * The SeqWeightTable class keeps track of sequence-specific bias parameters.
 * Allows for the bias associated with a given sequence to be calculated, and
 * for the bias parameters to be updated based on additional observations.  All
 * values are stored in log space.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class SeqWeightTable {
  /**
   * A private size_t specifying the order of the Markov chains used to model
   * the sequences.
   */
  size_t _order;
  /**
   * A private MarkovModel that stores the observed conditional nucleotide
   * frequencies (logged) in the bias window surrounding the fragment end.
   */
  MarkovModel _observed;
  /**
   * A private MarkovModel that stores the expected condtional nucleotide
   * frequencies (logged) based on an assumption of equal abundance of all
   * targets.
   */
  MarkovModel _expected;
public:
  /**
   * SeqWeightTable Constructor.
   * @param window_size an unsigned integer specifying the size of the bias
   *        window surrounding fragment ends.
   * @param order a size_t specifying the order to use for the Markov chains
   *        modelling the sequence.
   * @param alpha a double specifying the strength of the uniform prior
   *        (logged pseudo-counts for each parameter).
   */
  SeqWeightTable(size_t window_size, size_t order, double alpha);
  /**
   * A second constructor that loads the distribution from a parameter file.
   * Note that the values should not be modified after using this constructor.
   * @param window_size an unsigned integer specifying the size of the bias
   *        window surrounding fragment ends. Must match file.
   * @param order a size_t specifying the order to use for the Markov chains
   *        modelling the sequence. Must match file.
   * @param param_file_name a string specifying the path to the parameter file.
   * @param identifier a string specifying the header for these parameters in
   *        the file.
   */
  SeqWeightTable(size_t window_size, size_t order, std::string param_file_name,
                 std::string identifier);
  /**
   * A member function that overwrites the "observed" parameters with those from
   * another SeqWeightTable.
   * @param other another SeqWeightTable from which to copy the parameters.
   */
  void copy_observed(const SeqWeightTable& other);
  /**
   * A member function that overwrites the "expected" parameters with those from
   * another SeqWeightTable.
   * @param other another SeqWeightTable from which to copy the parameters.
   */
  void copy_expected(const SeqWeightTable& other);
  /**
   * A member function that increments the expected counts for a sliding window
   * through the given target sequence by some mass.
   * @param seq the target sequence.
   * @param mass the amount to increment by in the parameter table.
   * @param fl_cdf the fragment length CDF.
   */
  void increment_expected(const Sequence& seq, double mass,
                          const std::vector<double>& fl_cdf);
  /**
   * A member function that normalizes the expected counts and fills in the
   * lower-ordered marginals.
   */
  void normalize_expected();
  /**
   * A member function that increments the observed counts for the given
   * fragment sequence by some (logged) mass.
   * @param seq the target sequence (possibly reverse complemented) to which the
   *        fragment end maps.
   * @param i the index into the sequence at which to center the bias window, ie
   *        where the fragment starts/ends.
   * @param mass the amount to increment by (logged)
   */
   void increment_observed(const Sequence& seq, size_t i, double mass);
  /**
   * A member function that calculates the bias weight (logged) of a window.
   * This is the ratio of the observed and expected weights given by the two
   * Markov models.
   * @param seq the target sequence.
   * @param i the central point of the bias window, ie the fragment end.
   * @return The bias weight for the window.
   */
   double get_weight(const Sequence& seq, size_t i) const;
  /**
   * A member function that appends the marginal and conditional probabilities
   * for the foreground and background Markov models to the given file,
   * formatted in tables for easy readability.
   * @param outfile the file to append to.
   */
  void append_output(std::ofstream& outfile) const;
};

/**
 * The BiasBoss class keeps track of sequence-specific and positional bias.
 * It allows for the bias associated with a given fragment end to be calculated, and

 * for the bias parameters to be updated based on additional observations.  All stored
 * and returned values are in log space.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class BiasBoss
{
  /**
   * A private size_t specifying the order of the Markov chains used to model
   * the sequences.
   */
  size_t _order;
  /**
   * A private SeqWeightTable that stores the 5' sequence-specific bias
   * parameters (logged).
   */
  SeqWeightTable _5_seq_bias;
  /**
   * a private SeqWeightTable that stores the 3' sequence-specific bias
   * parameters (logged).
   */
  SeqWeightTable _3_seq_bias;

public:
  /**
   * BiasBoss Constructor.
   * @param order a size_t specifying the order of the Markov chains used to
   *        model the sequences.
   * @param alpha a double specifying the strength of the uniform prior (logged
   *        pseudo-counts for each parameter).
   */
  BiasBoss(size_t order, double alpha);
  /**
   * A second constructor that loads the distributions from a parameter file.
   * Note that the values should not be modified after using this constructor.
   * @param order a size_t specifying the order to use for the Markov chains
   *        modelling the sequence. Must match file.
   * @param param_file_name a string specifying the path to the parameter file.
   */  //TODO: Detect order from file.
  BiasBoss(size_t order, std::string param_file_name);
  /**
   * An accessor for the order of the Markov chains used to model the sequences.
   * @return The order of the Markov chains used to model the sequences.
   */
  size_t order() const { return _order; }
  /**
   * A member function that copies the observed parameters from another
   * BiasBoss.
   * @param other a BiasBoss to copy the parameters from.
   */
  void copy_observations(const BiasBoss& other);
  /**
   * A member function that copies the expected parameters from another
   * BiasBoss.
   * @param other a BiasBoss to copy the parameters from.
   */
  void copy_expectations(const BiasBoss& other);
  /**
   * A member function that updates the expectation parameters assuming uniform
   * abundance of and coverage accross the target's sequence.
   * @param targ the target to measure expected counts from
   */
   void update_expectations(const Target& targ,
                            double mass = 0,
                            const std::vector<double>& fl_cdf = std::vector<double>());
  /**
   * A member function that normalizes the expected counts and fills in the
   * lower-ordered marginals.
   */
  void normalize_expectations();
  /**
   * A member function that updates the observed parameters given a fragment
   * mapping to a target and its logged probabilistic assignment value.
   * @param hit the fragment hit (alignment).
   * @param mass the logged probabality of the mapping, which is the amount to
   *        increment the observed counts by.
   */
  void update_observed(const FragHit& hit, double mass);
  /**
   * A member function that returns the 5' and 3' bias values at each position
   * in a given target based on the current bias parameters.
   * @param start_bias a vector containing the logged bias for each 5' start
   *        site in the target.
   * @param end_bias a vector containing the logged bias for each 3' end site in
   *        the target.
   * @param targ the target for which to calculate the bias.
   * @return The product of the average 5' and 3' bias (logged).
   */
  double get_target_bias(std::vector<float>& start_bias,
                         std::vector<float>& end_bias,
                         const Target& targ) const;
  /**
   * A member function that appends the 5' and 3' bias parameters to the given
   * file, formatted in tables for easy readability.
   * @param outfile the file to append to.
   */
  void append_output(std::ofstream& outfile) const;
};

#endif

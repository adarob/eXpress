/**
 * fld.h
 *  express
 *
 *  Created by Adam Roberts on 3/20/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 */

#ifndef FLD_H
#define FLD_H

#include <vector>
#include <string>

/**
 * The FLD class keeps track of the observed fragment length distribution. It is
 * initialized with a Gaussian prior with parameters specified by the arguments
 * to the constructor.  A small "Gaussian" kernel is then added for each
 * observation. All mass values and probabilities are stored and returned in log
 * space (except in to_string).
 */
class FLD {
  /**
   * A private vector that stores the observed (logged) mass for each length.
   */
  std::vector<double> _hist;
  /**
   * A private double that stores the total observed (logged) mass.
   */
  double _tot_mass;
  /**
   * A private double that stores the (logged) sum of the product of observed
   * lengths and masses for quick mean calculations.
   */
  double _sum;
  /**
   * A private int that stores the minimum observed length.
   */
  size_t _min;
    
public:
  /**
   * FLD Constructor.
   * @param alpha double that sets the average pseudo-counts (logged).
   * @param max_val an integer that sets the maximum allowable length.
   * @param mean a size_t for the mean of the prior gaussian distribution.
   * @param std_dev a size_t for the standard deviation of the prior gaussian
   *        distribution.
   */
  FLD(double alpha, size_t max_val, size_t mean, size_t std_dev);
  /**
   * An accessor for the maximum allowed length.
   * @return Max allowed length.
   */
  size_t max_val() const;
  /**
   * An accessor for the minimum observed length (1 initially).
   * @return Minimum observed length.
   */
  size_t min_val() const;
  /**
   * An accessor for the mean length in the distribution.
   * @return Mean observed length.
   */
  double mean() const;
  /**
   * A member function that updates the distribution based on a new length
   * observation.
   * @param len an integer for the observed length.
   * @param mass a double for the mass (logged) to add.
   */
  void add_val(size_t len, double mass);
  /**
   * An accessor for the (logged) probability of a given length.
   * @param len an integer for the length to return the probability of.
   * @return (logged) probability of observing the given length.
   */
  double pmf(size_t len) const;
  /**
   * A member function that returns a vector containing the (logged) cumulative
   * mass function.
   * @return (Logged) cmf of fragment lengths.
   */
  std::vector<double> cmf() const;
  /**
   * An accessor for the (logged) observation mass (including pseudo-counts).
   * @return Total observation mass.
   */
  double tot_mass() const;
  /**
   * A member function that returns a string containing the current
   * distribution.
   * @return Space-separated string of probabilities ordered from length 0 to
   *         max_val (non-logged).
   */
  std::string to_string() const;
  /**
   * A member function that appends the FLD parameters to the end of the given
   * file.
   * @param outfile the file to append to.
   */
  void append_output(std::ofstream& outfile) const;
};

#endif

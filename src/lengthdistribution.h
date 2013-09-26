/**
 *  lengthdistribution.h
 *  express
 *
 *  Created by Adam Roberts on 1/30/13.
 *  Copyright 2013 Adam Roberts. All rights reserved.
 */

#ifndef LengthDistribution_H
#define LengthDistribution_H

#include <vector>
#include <string>

/**
 * The LengthDistribution class keeps track of the observed length distribution.
 * It is initialized with a Gaussian prior with parameters specified by the
 * arguments to the constructor. An argument-specified binomial kernel is then
 * added for each observation. All mass values and probabilities are stored and
 * returned in log space (except in to_string).
 */
class LengthDistribution {
  /**
   * A private vector that stores the (logged) kernel values.
   **/
  std::vector<double> _kernel;
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
  /**
   * A size for internal binning of the lengths in the distribution.
   */
  size_t _bin_size;
  
public:
  /**
   * LengthDistribution Constructor.
   * @param alpha double that sets the average pseudo-counts (logged).
   * @param max_val an integer that sets the maximum allowable length.
   * @param prior_mu a size_t for the mean of the prior gaussian distribution.
            If 0, a uniform distribution is used instead.
   * @param prior_sigma a size_t for the standard deviation of the prior
   *        gaussian distribution.
   * @param kernel_n a size_t specifying the number of trials in the kernel
   *        binomial distribution. Must be odd.
   * @param kernel_p a double specifying the success probability for the kernel
   *        binomial distribution.
   * @param bin_size a size_t specifying the size of bins to use internally to
   *        reduce the number of parameters in the distribution.
   */
  LengthDistribution(double alpha, size_t max_val, size_t prior_mu,
                     size_t prior_sigma, size_t kernel_n, double kernel_p,
                     size_t bin_size = 1);
  /**
   * A second constructor that loads the distribution from a parameter file.
   * Note that the values should not be modified (add_val should not be called)
   * after using this constructor. The bin size is set to 1.
   * @param param_file_name a string specifying the path to the parameter file.
   * @param length_type a string specifying the type of length distribution
   *        to be matched in the parameter file.
   */
  LengthDistribution(std::string param_file_name, std::string length_type);
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
   * A member function that returns a (logged) cumulative mass for a given
   * length.
   * @param len an integer for the length to return the cmf value of.
   * @return (Logged) cmf value of length.
   */
  double cmf(size_t len) const;
  /**
   * A member function that returns a vector containing the (logged) cumulative
   * mass function *for the bins*.
   * @return (Logged) cmf of bins.
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
   * A member function that appends the LengthDistribution parameters to the end
   * of the given file.
   * @param outfile the file to append to.
   * @param length_type a string specifying the type of length the distribution
   *        is of (ie. "Fragment" or "Target") to be included in the header.
   */
  void append_output(std::ofstream& outfile, std::string length_type) const;
};

#endif

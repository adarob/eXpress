#ifndef BIASCORRECTION_H
#define BIASCORRECTION_H

/*
 *  biascorrection.h
 *  express
 *
 *  Created by Adam Roberts on 4/5/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 *
 */

#include <vector>
#include <string>
#include "markovmodel.h"
#include "frequencymatrix.h"
#include "fld.h"

class Target;
class FragHit;

/**
 * The SeqWeightTable class keeps track of sequence-specific bias parameters.
 * It allows for the bias associated with a given sequence to be calculated, and for the bias
 * parameters to be updated based on additional observations.  All values stored in log space.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class SeqWeightTable
{
    /**
     * a private MarkovModel that stores the observed conditional nucleotide frequencies (logged) in the bias window
     * surrounding the fragment end
     */
    MarkovModel _observed;
    
    /**
     * a private MarkovModel that stores the expected condtional nucleotide frequencies (logged) based on an assumption
     * of equal expression of all targets
     */
    MarkovModel _expected;
    
public:
    
    /**
     * SeqWeightTable Constructor
     * @param window_size an unsigned integer specifying the size of the bias window surrounding fragment ends
     * @param alpha a double specifying the strength of the uniform prior (logged pseudo-counts for each parameter)
     */
    SeqWeightTable(size_t window_size, double alpha);
    
    
    /**
     * a member function that overwrites the "observed" parameters with those from another SeqWeightTable
     * @param other another SeqWeightTable from which to copy the parameters
     */
    void copy_observed(const SeqWeightTable& other);

    /**
     * a member function that overwrites the "expected" parameters with those from another SeqWeightTable
     * @param other another SeqWeightTable from which to copy the parameters
     */    
    void copy_expected(const SeqWeightTable& other);
    
    /**
     * a member function that increments the expected counts for a sliding window through
     * the given target sequence by some mass
     * @param seq the target sequence
     * @param mass the amount of used to weight the target's sequence in the parameter table
     */
    void increment_expected(const Sequence& seq, double mass, const std::vector<double>& fl_cdf); 
    
    /**
     * a member function that normalizes the expected counts and converts them to the log scale
     */
    void normalize_expected();
    
    /**
     * a member function that increments the observed counts for the given fragment sequence by some mass (logged)
     * @param seq the target sequence (possibly reverse complemented) to which the fragment end maps
     * @param i the index into the sequence at which to center the bias window (where the fragment starts/ends)
     * @param mass the fragment's mass
     */
    void increment_observed(const Sequence& seq, size_t i, double mass);
    
    /**
     * a member function that calculates the bias weight (logged) of a bias window
     * @param seq the target sequence the fragment hits to
     * @param i the fragment end point (the central point of the bias window)
     * @return the bias weight for the bias window which is the product of the individual nucleotide bias weights
     */
    double get_weight(const Sequence& seq, size_t i) const;
    
    /**
     * a member function that returns a string containing the positional nucleotide probabilities in column-major order (A,C,G,T)
     * @return the string representation of the positional nucleotide probabilities
     */
    std::string to_string() const;
    
    /**
     * a member function that outputs the positional nucleotide probabilities in matrix format with nucleotides (A,C,G,T) as
     * rows and window position as columns
     * @param outfile the file to append to
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
     * a private SeqWeightTable that stores the 5' sequence-specific bias parameters (logged)
     */
    SeqWeightTable _5_seq_bias;
    
    /**
     * a private SeqWeightTable that stores the 3' sequence-specific bias parameters (logged)
     */
    SeqWeightTable _3_seq_bias;
    
  public:
    
    /**
     * BiasBoss Constructor
     * @param alpha a double specifying the strength of the uniform prior (logged pseudo-counts for each parameter)
     */
    BiasBoss(double alpha);
        
    //DOC
    void copy_observations(const BiasBoss& other);
    void copy_expectations(const BiasBoss& other);
    
    /**
     * a member function that updates the expectation parameters (sequence-specific and positional)
     * assuming uniform expression of and accross the target's sequence
     * @param targ the target to measure expected counts from
     */
    void update_expectations(const Target& targ, double mass=0, const std::vector<double>& fl_cdf=std::vector<double>());
    
    /**
     * a member function that normalizes the expected counts and converts them to the log scale
     */
    void normalize_expectations();
    
    /**
     * a member function that updates the observed parameters (sequence-specific and positional) 
     * given a fragment mapping to a target and its logged probabilistic assignment
     * @param hit the fragment hit (alignment)
     * @param mass the logged probabality of the mapping, which is the amount to update the observed counts by
     */
    void update_observed(const FragHit& hit, double mass);
    
    /**
     * a member function that returns the 5' and 3' bias values at each position in a given target
     * based on the current bias parameters
     * @param start_bias a vector containing the logged bias for each 5' start site in the target
     * @param end_bias a vector containing the logged bias for each 3' end site in the target
     * @param targ the target for which to calculate the logged bias
     * @return the product of the average 5' and 3' bias (logged)
     */
    double get_target_bias(std::vector<float>& start_bias, std::vector<float>& end_bias, const Target& targ) const;
    
    /**
     * a member function that returns a string containing the observed positional nucleotide probabilities (non-logged) in column-major order (A,C,G,T)
     * @return the string representation of the observed probabilities
     */
    std::string to_string() const;
    
    /**
     * a member function that outputs the positional and sequence-specific bias parameter matrices
     * @param outfile the file to append to
     */
    void append_output(std::ofstream& outfile) const;
};


#endif

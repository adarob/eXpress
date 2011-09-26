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
#include <boost/thread.hpp>
#include "frequencymatrix.h"

class Transcript;
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
     * a private FrequencyMatrix that stores the observed nucleotide frequencies (logged) in the bias window
     * surrounding the fragment end
     */
    FrequencyMatrix _observed;
    
    /**
     * a private FrequencyMatrix that stores the expected nucleotide frequencies (logged) based on an assumption
     * of equal expression of all transcripts
     */
    FrequencyMatrix _expected;
    
    /**
     * a private mutex block access for multi-threaded bias updating
     */
    mutable boost::mutex _lock;
    
public:
    
    /**
     * SeqWeightTable Constructor
     * @param window_size an unsigned integer specifying the size of the bias window surrounding fragment ends
     * @param alpha a double specifying the strength of the uniform prior (logged pseudo-counts for each paramater)
     */
    SeqWeightTable(size_t window_size, double alpha);
    
    /**
     * a member function that increments the expected counts for the given nucleotide by 1 (logged)
     * @param c a char representing a nucleotide that has been observed in the transcriptome
     */
    void increment_expected(char c); 
    
    /**
     * a member function that increments the observed counts for the given fragment sequence by some mass (logged)
     * @param seq a string of nucleotides in the bias window for the sequenced fragment end
     * @param normalized_mass the mass (logged probabilistic assignment) of the fragment normalized by its estimated expression
     */
    void increment_observed(std::string& seq, double normalized_mass);
    
    /**
     * a member function that calculates the bias weight (logged) of a bias window
     * @param seq the transcript sequence the fragment hits to
     * @param i the fragment end point (the central point of the bias window)
     * @return the bias weight for the bias window which is the product of the individual nucleotide bias weights
     */
    double get_weight(const std::string& seq, int i) const;
    
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
 * The PosWeightTable class keeps track of fractional position bias parameters in log space.
 * It allows for the bias associated with a given fractional position to be calculated, and for the bias
 * parameters to be updated based on additional fragment observations.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class PosWeightTable
{
    /**
     * a private FrequencyMatrix that stores the fragment frequencies (logged) in different positional bins
     */
    FrequencyMatrix _observed;
    
    /**
     * a private FrequencyMatrix that stores the expected fragment frequencies (logged) in different positional
     * bins based on an assumption of equal expression of all transcripts
     */
    FrequencyMatrix _expected;
    
    /**
     * a private vector of unsigned integers specifying the bin ranges for transcript lengths
     */
    const std::vector<size_t> _len_bins;
    
    /**
     * a private vector of doubles specifying the bin ranges for fractional positions
     */
    const std::vector<double> _pos_bins;
    
    /**
     * a private mutex block access for multi-threaded bias updating
     */
    mutable boost::mutex _lock;
    
public:
    
    /**
     * PosWeightTable Constructor
     * @param len_bins a vector of unsigned integers specifying the bin ranges for transcript lengths
     * @param pos_bins a vector of doubles specifying the bin ranges for fractional positions
     * @param alpha a double specifying the strength of the uniform prior (logged pseudo-counts for each paramater)
     */
    PosWeightTable(const std::vector<size_t>& len_bins, const std::vector<double>& pos_bins, double alpha);
    
    const std::vector<size_t>& len_bins() const { return _len_bins; }
    const std::vector<double>& pos_bins() const { return _pos_bins; }

    /**
     * a member function that increments the expected counts for the given fractional position by 1 (logged)
     * @param len the transcript length
     * @param pos the fractional transcript position
     */
    void increment_expected(size_t len, double pos); 
    
    /**
     * a member function that increments the expected counts for the given fractional position bin by 1 (logged)
     * @param l the transcript length bin
     * @param p the fractional transcript position bin
     */
    void increment_expected(size_t l, size_t p); 

    
    /**
     * a member function that increments the observed counts for the given fragment position by some mass (logged)
     * @param len the transcript length
     * @param pos the fractional transcript position
     * @param normalized_mass the mass (logged probabilistic assignment) of the fragment normalized by its estimated expression
     */
    void increment_observed(size_t len, double pos, double normalized_mass);
    
    /**
     * a member function that increments the observed counts for the given fragment position bin by some mass (logged)
     * @param l the transcript length bin
     * @param p the fractional transcript position bin
     * @param normalized_mass the mass (logged probabilistic assignment) of the fragment normalized by its estimated expression
     */
    void increment_observed(size_t l, size_t p, double normalized_mass);

    /**
     * a member function that return the bias weight (logged) of a fractional transcript position
     * @param len the transcript length
     * @param pos the fractional transcript position
     * @return the logged bias weight for the fractional transcript position
     */
    double get_weight(size_t len, double pos) const;
    
    /**
     * a member function that return the bias weight (logged) of a fractional transcript position bin
     * @param l the transcript length bin
     * @param p the fractional transcript position bin
     * @return the logged bias weight for the fractional transcript position
     */
    double get_weight(size_t l, size_t p) const;

    
    /**
     * a member function that outputs the fractional position probabilities in matrix format with length bins
     * as rows and fractional position bins as columns
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
    
    /**
     * a private SeqWeightTable that stores the 5' sequence-specific bias parameters (logged)
     */
    PosWeightTable _5_pos_bias;
    
    /**
     * a private SeqWeightTable that stores the 3' sequence-specific bias parameters (logged)
     */
    PosWeightTable _3_pos_bias;
    
public:
    
    /**
     * BiasBoss Constructor
     * @param alpha a double specifying the strength of the uniform prior (logged pseudo-counts for each paramater)
     */
    BiasBoss(double alpha);
    
    /**
     * a member function that updates the expectation parameters (sequence-specific and positional)
     * assuming uniform expression of and accross the transcript's sequence
     * @param trans the transcript to measure expected counts from
     */
    void update_expectations(const Transcript& trans);
    
    /**
     * a member function that updates the observed parameters (sequence-specific and positional) 
     * given a fragment mapping to a transcript and its logged probabilistic assignment
     * @param hit the fragment hit (alignment)
     * @param mass the logged probabality of the mapping, which is the amount to update the observed counts by
     */
    void update_observed(const FragHit& hit, double mass);
    
    /**
     * a member function that returns the 5' and 3' bias values at each position in a given transcript
     * based on the current bias parameters
     * @param start_bias a vector containing the logged bias for each 5' start site in the transcript
     * @param end_bias a vector containing the logged bias for each 3' end site in the transcript
     * @param trans the transcript for which to calculate the logged bias
     * @return the product of the average 5' and 3' bias (logged)
     */
    double get_transcript_bias(std::vector<double>& start_bias, std::vector<double>& end_bias, const Transcript& trans) const;
    
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

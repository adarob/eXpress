#ifndef BIASCORRECTION_H
#define BIASCORRECTION_H

/*
 *  biascorrection.h
 *  cufflinks
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
class FragMap;

/**
 * SeqWeightTable class.  This class keeps track of sequence-specific bias parameters.
 * It allows for the bias associated with a given sequence to be calculated, and for the bias
 * parameters to be updated based on additional observations.
 */
class SeqWeightTable
{
    /**
     * a private FrequencyMatrix that stores the observed nucleotide frequencies in the bias window
     * surrounding the fragment end
     */
    FrequencyMatrix _observed;
    
    /**
     * a private FrequencyMatrix that stores the expected nucleotide frequencies based on an assumption
     * of equal expression of all transcripts
     */
    FrequencyMatrix _expected;
    
    /**
     * a private mutex block access for multi-threaded bias updating
     */
    mutable boost::mutex _lock;
    
    /**
     * a private function that converts a nucleotide character into an integer representation
     * @param c a char representing a nucleotide
     * @return an unsigned integer representing the nucleotide (0-3) or 4 if not a nucleotide
     */
    size_t ctoi(char c) const;
public:
    
    /**
     * SeqWeightTable Constructor
     * @param window_size an unsigned integer specifying the size of the bias window surrounding fragment ends
     * @param alpha a double specifying the strength of the uniform prior (pseudo-counts for each paramater)
     */
    SeqWeightTable(size_t window_size, double alpha);
    
    /**
     * a member function that increments the expected counts for the given nucleotide by 1
     * @param c a char representing a nucleotide that has been observed in the transcriptome
     */
    void increment_expected(char c); 
    
    /**
     * a member function that increments the observed counts for the given fragment sequence by some mass
     * @param seq a string of nucleotides in the bias window for the sequenced fragment end
     * @param normalized_mass the mass (probabilistic assignment) of the fragment normalized by its estimated expression
     */
    void increment_observed(std::string& seq, double normalized_mass);
    
    /**
     * a member function that calculates the bias weight of a bias window
     * @param seq the transcript sequence the fragment maps to
     * @param i the fragment end point (the central point of the bias window)
     * @return the bias weight for the bias window which is the product of the individual nucleotide bias weights
     */
    double get_weight(const std::string& seq, size_t i) const;
};

/**
 * PosWeightTable class.  This class keeps track of fractional position bias parameters.
 * It allows for the bias associated with a given fractional position to be calculated, and for the bias
 * parameters to be updated based on additional fragment observations.
 */
class PosWeightTable
{
    /**
     * a private FrequencyMatrix that stores the fragment frequencies in different positional bins
     */
    FrequencyMatrix _observed;
    
    /**
     * a private FrequencyMatrix that stores the expected fragment frequencies in different positional
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
     * @param alpha a double specifying the strength of the uniform prior (pseudo-counts for each paramater)
     */
    PosWeightTable(const std::vector<size_t>& len_bins, const std::vector<double>& pos_bins, double alpha);
    
    const std::vector<size_t>& len_bins() const { return _len_bins; }
    const std::vector<double>& pos_bins() const { return _pos_bins; }

    /**
     * a member function that increments the expected counts for the given fractional position by 1
     * @param len the transcript length
     * @param pos the fractional transcript position
     */
    void increment_expected(size_t len, double pos); 
    
    /**
     * a member function that increments the expected counts for the given fractional position by 1
     * @param l the transcript length bin
     * @param p the fractional transcript position bin
     */
    void increment_expected(size_t l, size_t p); 

    
    /**
     * a member function that increments the observed counts for the given fragment position by some mass
     * @param len the transcript length
     * @param pos the fractional transcript position
     * @param normalized_mass the mass (probabilistic assignment) of the fragment normalized by its estimated expression
     */
    void increment_observed(size_t len, double pos, double normalized_mass);
    
    /**
     * a member function that increments the observed counts for the given fragment position by some mass
     * @param l the transcript length bin
     * @param p the fractional transcript position bin
     * @param normalized_mass the mass (probabilistic assignment) of the fragment normalized by its estimated expression
     */
    void increment_observed(size_t l, size_t p, double normalized_mass);

    /**
     * a member function that return the bias weight of a fractional transcript position
     * @param len the transcript length
     * @param pos the fractional transcript position
     * @return the bias weight for the fractional transcript position
     */
    double get_weight(size_t len, double pos) const;
    
    /**
     * a member function that return the bias weight of a fractional transcript position
     * @param l the transcript length bin
     * @param p the fractional transcript position bin
     * @return the bias weight for the fractional transcript position
     */
    double get_weight(size_t l, size_t p) const;

};

/**
 * BiasBoss class.  This class keeps track of sequence-specific and positional bias.
 * It allows for the bias associated with a given fragment end to be calculated, and 
 * for the bias parameters to be updated based on additional observations.
 */
class BiasBoss
{
    /**
     * a private SeqWeightTable that stores the 5' sequence-specific bias parameters
     */
    SeqWeightTable _5_seq_bias;
    
    /**
     * a private SeqWeightTable that stores the 3' sequence-specific bias parameters
     */
    SeqWeightTable _3_seq_bias;
    
    /**
     * a private SeqWeightTable that stores the 5' sequence-specific bias parameters
     */
    PosWeightTable _5_pos_bias;
    
    /**
     * a private SeqWeightTable that stores the 3' sequence-specific bias parameters
     */
    PosWeightTable _3_pos_bias;
    
public:
    
    /**
     * BiasBoss Constructor
     * @param alpha a double specifying the strength of the uniform prior (pseudo-counts for each paramater)
     */
    BiasBoss(double alpha);
    
    /**
     * a member function that updates the expectation parametersn(sequence-specific and positional)
     * assuming uniform expression of and accross the transcript's sequence
     * @param trans the transcript to measure expected counts from
     */
    void update_expectations(const Transcript& trans);
    
    /**
     * a member function that updates the observed parameters (sequence-specific and positional) 
     * given a fragment mapping to a transcript and its probabilistic assignment
     * @param frag the fragment mapping
     * @param trans the transcript mapped to by the fragment
     * @param mass the probabality of the mapping, which is the amount to update the observed counts by
     */
    void update_observed(const FragMap& frag, const Transcript& trans, double mass);
    
    /**
     * a member function that returns the 5' and 3' bias values at each position in a given transcript
     * based on the current bias parameters
     * @param start_bias a vector containing the bias for each 5' start site in the transcript
     * @param start_bias a vector containing the bias for each 3' end site in the transcript
     * @param trans the transcript for which to calculate the bias
     * @return the product of the average 5' and 3' bias
     */
    double get_transcript_bias(std::vector<double>& start_bias, std::vector<double>& end_bias, const Transcript& trans) const;
};


#endif

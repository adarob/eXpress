//
//  transcripts.h
//  express
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#ifndef TRANSCRIPTS_H
#define TRANSCRIPTS_H

#include <string>
#include <map>
#include <boost/unordered_map.hpp>
#include <vector>
#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include "main.h"
#include "bundles.h"

class FLD;
class FragHit;
class BiasBoss;
class MismatchTable;

typedef size_t TransID;


/**
 * The Transcript class is used to store objects for the transcripts
 * being mapped to.  Besides storing basic information about the object
 * (id, length), it also stores a mass based on the number of 
 * fragments mapping to the object.  To help with updating this number,
 * it returns the likelihood that a given fragment originated from it.
 * These values are stored and returned in log space.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class Transcript
{
    /**
     * a private pointer to the struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
     */
    const Globals* _globs;
    
    /**
     * a private TransID that stores the hashed transcript name
     */
    TransID _id;
    
    /**
     * a private string that stores the transcript name
     */
    std::string _name;
    
    /**
     * a private string that stores the transcript sequence
     */
    std::string _seq;
    
    double _alpha;
    
    /**
     * a private double that stores the (logged) mass based on observed fragment mapping probabilities
     */
    double _mass;
    
    //FIX
    double _ambig_mass;

    /**
     * a private double that stores the (logged) binomal variance of the mass
     */
    double _binom_var;
    
    //FIX
    double _samp_var;
    
    /**
     * a private double that stores the total (logged) mass of all fragments mapped to the transcript
     */
    double _tot_mass;
    
    //FIX
    double _tot_ambig_mass;
    
    //FIX
    double _tot_unc;
    
    /**
     * a private double that stores the estimated counts of observed fragment mappings
     * used after an initial run of the online EM to improve estimates
     */
    double _est_counts;
    
    /**
     * a private double that stores the variance of the estimated counts
     * used after an initial run of the online EM to improve estimates
     */
    double _est_counts_var;
    
    /**
     * a private size_t that stores the number of fragments (non-logged) uniquely mapping
     * to this transcript
     */
    size_t _uniq_counts;
    
    /**
     * a private size_t that stores the fragment counts (non-logged) for the bundle
     * the total bundle counts is the sum of this value for all transcripts in the bundle
     */
    size_t _tot_counts;
    
    /**
     * a private pointer to the Bundle this Transcript is a member of
     */
    Bundle* _bundle;
    
    /**
     * a private mutex to provide thread-safety for bias variables with threaded update
     */
    mutable boost::mutex _bias_lock;
    
    /**
     * a point to a private double vector storing the (logged) 5' bias at each position (accessed by multiple threads)
     */
    std::vector<float>* _start_bias;
    
    /**
     * a pointer to a private double vector storing the (logged) 3' bias at each position (accessed by multiple threads)
     */
    std::vector<float>* _end_bias;
    
    /**
     * a private double storing the (logged) product of the average 3' and 5' biases for the transcript (accessed by multiple threads)
     */
    double _avg_bias;
    
    /**
     * a private double storing the most recently updated (logged) effective length as calculated by the bias updater thread
     */
    double _cached_eff_len;
    
public:
    
    /**
     * Transcript Constructor
     * @param name a string that stores the transcript name
     * @param seq a string that stores the transcript sequence
     * @param alpha a double that specifies the intial pseudo-counts (non-logged)
     * @param globs a pointer to the struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
     */
    Transcript(const size_t id, const std::string& name, const std::string& seq, double alpha, const Globals* globs);
    
    ~Transcript()
    {
        if (_start_bias)
            delete _start_bias;
        if (_end_bias)
            delete _end_bias;
    }
    /**
     * a member function that returns the transcript name
     * @return string containing transcript name
     */
    const std::string& name() const { return _name; }
    
    /**
     * a member function that returns the transcript id
     * @return TransID transcript ID
     */
    TransID id() const { return _id; }
    /**
     * a member function that returns the transcript sequence
     * @return string containing transcript sequence
     */
    const std::string& seq() const { return _seq; }
    
    /**
     * a member function that returns the transcript length
     * @return transcript length
     */
    size_t length() const { return _seq.length(); }
    
    /**
     * a member function that returns the current (logged) probabilistically assigned fragment mass
     * FIX
     * @return logged mass
     */
    double mass(bool with_pseudo = false) const;
    
    //FIX
    double ambig_mass() const { return _ambig_mass; }
    
    /**
     * a member function that returns the total (logged) mass of all fragments mapped to the transcript
     * @return logged total mass
     */
    double tot_mass() const { return _tot_mass; }
    
    //FIX
    double tot_ambig_mass() const { return _tot_ambig_mass; }
    
    /**
     * a member function that returns the current (logged) binomial variance
     * @return logged mass variance
     */
    double binom_var() const { return _binom_var; }
    
    //FIX
    double samp_var() const { return _samp_var; }
    
    //FIX
    double tot_uncertainty() const { return _tot_unc; }
    
    /**
     * a member function that returns the current (logged) estimated counts
     * not valid in the first online EM round
     * @return estimated counts
     */
    double est_counts() const { return _est_counts; }
         
    /**
     * a member function that returns the current (logged) variance on estimated counts
     * not valid in the first online EM round
     * @return variance estimated counts
     */
    double est_counts_var() const { return _est_counts_var; }

    /**
     * a member function that readies the transcript object for the next round of batch EM
     */
    void round_reset();
    
    /**
     * a member function that returns the current count of fragments mapped to this transcript (uniquely or ambiguously)
     * @return total fragment count
     */
    size_t tot_counts() const { return _tot_counts; }
    
    /**
     * a member function that returns the current count of fragments uniquely mapped to this transcript
     * @return unique fragment count
     */
    size_t uniq_counts() const { return _uniq_counts; }
    
    /**
     * a member function that returns the Bundle this Transcript is a member of
     * @return a pointer to the Bundle this transcript is a member of
     */
    Bundle* bundle() { return _bundle; }

    /**
     * a member function that set the Bundle this Transcript is a member of
     * @param b a pointer to the Bundle to set this Transcript as a member of
     */
    void bundle(Bundle* b) { _bundle = b; }
    
    /**
     * a member function that increases the expected fragment counts and variance by a given (logged) fragment mass
     * @param p a double for the (logged) probability that the fragment was generated by this transcript
     * @param mass a double specifying the (logged) mass of the fragment being mapped
     */
    //FIX
    void add_mass(double p, double v, double mass);
    
    /**
     * a member function that increases the estimated counts and estimated count variance based on the probabilistic
     * assignment of a fragment
     * @param p a double for the (logged) probability that the fragment was generated by this transcript
     */
    void add_prob_count(double p);

    /**
     * a member function that increases the count of fragments mapped to this transcript
     * @param uniq a bool specifying whether or not the fragment uniquely maps to this transcript
     * @param incr_amt a size_t to increase the counts by
     */
    void incr_counts(bool uniq, size_t incr_amt = 1)
    {
        _tot_counts += incr_amt;
        _uniq_counts += incr_amt*uniq;
    }
    
    /**
     * a member function that returns (a value proportional to) the log likelihood the given fragment
     * originated from this transcript
     * @param frag a FragHit to return the likelihood of being originated from this transcript
     * @param with_pseudo a FragHit specifying whether or not pseudo-counts should be included in the calculation
     * @return (a value proportional to) the log likelihood the given fragment originated from this transcript
     */
    double log_likelihood(const FragHit& frag, bool with_pseudo) const;
    
    /**
     * a member function that calcualtes and returns the effective length of the transcript (logged)
     * @return the effective length of the transcript calculated as \f$ \tilde{l} = \sum_{l=1}^{L(t)}\sum_{i=1}^{L(t)} D(l)b_5[i]*b_3[i+l] \f$
     */
    //double effective_length() const;
    
    /**
     * a member function that calcualtes and returns the estimated effective length of the transcript (logged) using the avg bias
     * @return the estimated effective length of the transcript calculated as \f$ \tilde{l} = \bar{bias}\sum_{l=1}^{L(t)} D(l)(L(t) - l + 1) \f$
     */
    double est_effective_length() const;
    
    /**
     * a member function that returns the most recently estimated effective length (logged) as calculated by the bias updater thread 
     * @return the cached effective length of the transcript calculated
     */
    double cached_effective_length() const;

    /**
     * a member function that causes the transcript bias to be re-calculated by the _bias_table based on curent parameters
     */
    void update_transcript_bias();
    
};

typedef std::vector<Transcript*> TransMap;
typedef boost::unordered_map<std::string, size_t> TransIndex;
typedef boost::unordered_map<size_t, double> CovarMap;

/**
 * The TranscriptTable class is used to keep track of the Transcript objects for a run.
 * The constructor parses a fasta file to generate the Transcript objects and store them in a map
 * that allows them to be looked up based on their string id.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class TranscriptTable
{
    /**
     * a private pointer to the struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
     */
    const Globals* _globs;
        
    /**
     * a private map to look up pointers to Transcript objects by their TransID id
     */
    TransMap _trans_map;
    
    /**
     * the private table to keep track of Bundles
     */
    BundleTable _bundle_table;
    
    /**
     * a private CovarTable to look up the covariance for pairs of Transcripts by their combined hashed TransIDs
     * these values are stored positive, even though they are negative
     */
    CovarTable _covar_table;
    
    /**
     * a private double specifying the initial Transcript mass pseudo-counts for each bp (non-logged)
     */
    double _alpha;
    
    /**
     * a private function that validates and adds a transcript pointer to the table
     * @param name the name of the trancript
     * @param seq the sequence of the transcript
     * @param trans_index the transcript-to-index map from the alignment file
     * @param trans_lengths the transcript-to-length map from the alignment file (for validation)
     */
    void add_trans(const std::string& name, const std::string& seq, const TransIndex& trans_index, const TransIndex& trans_lengths);  
    
public:
    /**
     * TranscriptTable Constructor
     * @param trans_fasta_file a string storing the path to the fasta file from which to load transcripts
     * @param trans_index the transcript-to-index map from the alignment file
     * @param trans_lengths the transcript-to-length map from the alignment file
     * @param alpha a double that specifies the intial pseudo-counts for each bp of the transcripts (non-logged)
     * @param globs a pointer to the struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
     */
    TranscriptTable(const std::string& trans_fasta_file, const TransIndex& trans_index, const TransIndex& trans_lengths, double alpha, const Globals* globs);
    
    /**
     * TranscriptTable Destructor
     * deletes all of the transcript objects in the table
     */
    ~TranscriptTable();
    
    /**
     * a member function that returns a pointer to the transcript with the given id 
     * @param id of the transcript queried
     * @return pointer to the transcript wit the given id
     */
    Transcript* get_trans(TransID id);
    
    /**
     * a member function that readies all Transcript objects in the table for the next round of batch EM
     */
    void round_reset();
    
    /**
     * a member function that returns the number of transcripts in the table
     * @return number of transcripts in the table
     */
    size_t size() const { return _trans_map.size(); }
    
    /**
     * a member function that increases the covariance between two transcripts by the specified amount
     * these values are stored positive even though they are negative (logged)
     * @param trans1 one of the transcripts in the pair
     * @param trans2 the other transcript in the pair
     * @param covar a double specifying the amount to increase the pair's covariance by (logged)
     */
    void update_covar(TransID trans1, TransID trans2, double covar) { _covar_table.increment(trans1, trans2, covar); }
    
    /**
     * a member function that returns the covariance between two transcripts
     * these returned value will be the log of the negative of the true value
     * @param trans1 one of the transcripts in the pair
     * @param trans2 the other transcript in the pair
     * @return the negative of the pair's covariance (logged)
     */
    double get_covar(TransID trans1, TransID trans2) { return _covar_table.get(trans1, trans2); }
    
    /**
     * a member function that returns the number of pairs of transcripts with non-zero covariance
     * @return the number of transcript pairs with non-zero covariance
     */
    size_t covar_size() const { return _covar_table.size(); }
       
    /**
     * a member function that merges the given Bundles
     * @param b1 a pointer to the first Bundle to merge
     * @param b2 a pointer to the second Bundle to merge
     * @return a pointer to the merged Bundle
     */
    Bundle* merge_bundles(Bundle* b1, Bundle* b2);
    
    /**
     * a member function that returns the number of bundles in the partition
     * @return the number of bundles in the partition
     */
    size_t num_bundles();
    
    /**
     * a member function for driving a thread that continuously updates the transcript bias values
     */
    void threaded_bias_update();
    
    /**
     * a member function that outputs the final expression data in a file called 'results.xprs'
     * and (optionally) the variance-covariance matrix in 'varcov.xprs' in the given output directory
     * @param output_dir the directory to output the expression file to
     * @param tot_counts the total number of observed mapped fragments
     * @param output_varcov boolean specifying whether to also output the variance-covariance matrix
     */
    void output_results(std::string output_dir, size_t tot_counts, bool output_varcov);
};

#endif
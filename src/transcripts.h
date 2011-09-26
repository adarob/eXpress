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
#include <boost/functional/hash.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <vector>
#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include "main.h"


class FLD;
class FragHit;
class BiasBoss;
class MismatchTable;

typedef size_t TransID;

/**
 * a global function to hash transcript names to TransIDs
 */
static boost::hash<std::string> hash_trans_name;

/**
 * a global function to combine TransIDs into paired hashes
 */
static size_t hash_trans_pair(TransID trans1, TransID trans2)
{
    size_t seed = std::min(trans1, trans2);
    boost::hash_combine(seed, std::max(trans1, trans2));
    return seed;
}

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
     * a private string that stores the transcript name
     */
    std::string _name;
    
    /**
     * a private TransID that stores the hashed transcript name
     */
    TransID _id;
    
    /**
     * a private string that stores the transcript sequence
     */
    std::string _seq;
    
    /**
     * a private double that stores the (logged) mass based on observed fragment mappings
     */
    double _mass;

    /**
     * a private double that stores the (logged) variance of the mass
     */
    double _var;
    
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
     * a private size_t that stores a portion of the fragment counts (non-logged) for the bundle
     * the total bundle counts is the sum of this value for all transcripts in the bundle
     */
    size_t _bundle_counts;
    
    /**
     * a private mutex to provide thread-safety for bias variables with threaded update
     */
    mutable boost::mutex _bias_lock;
    
    /**
     * a private double vector storing the (logged) 5' bias at each position (accessed by multiple threads)
     */
    std::vector<double> _start_bias;
    
    /**
     * a private double vector storing the (logged) 3' bias at each position (accessed by multiple threads)
     */
    std::vector<double> _end_bias;
    
    /**
     * a private double storing the (logged) product of the average 3' and 5' biases for the transcript (accessed by multiple threads)
     */
    double _avg_bias;
    
    /**
     * a private double storing the (non-logged) initial effective length of the transcript (based on the prior FLD), ignoring bias
     */
    double _ub_eff_len;
    
    /**
     * a private pointer to the "global" Fragment Length Distribution (FLD) object
     */
    const FLD* _fld;
    
    /**
     * a private pointer to the "global" bias table object
     */
    const BiasBoss* _bias_table;
    
    /**
     * a private pointer to the "global" mismatch (error) table object
     */
    const MismatchTable* _mismatch_table;
    
    /**
     * a private member function that returns the (logged) total bias for a given fragment length
     * @param l a size_t that specifies the fragment length
     * @return the total bias for length (\f$ B = \bar{bias}(L(t) - l + 1) \f$)
     */
    double total_bias_for_length(size_t l) const;
    
public:
    
    /**
     * Transcript Constructor
     * @param name a string that stores the transcript name
     * @param seq a string that stores the transcript sequence
     * @param alpha a double that specifies the intial pseudo-counts (non-logged)
     * @param fld a pointer to the global Fragment Length Distribution (FLD) object
     * @param bias_table a pointer to the global BiasBoss object
     * @param mismatch_table a pointer to the global MismatchTable object
     */
    Transcript(const std::string& name, const std::string& seq, double alpha, const FLD* fld, const BiasBoss* bias_table, const MismatchTable* mismatch_table);
    
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
     * a member function that returns the current (logged) fragment mass
     * @return logged mass
     */
    double mass() const { return _mass; }
    
    /**
     * a member function that returns the current (logged) variance
     * @return logged mass variance
     */
    double var() const { return _var; }
    
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
     * a member function that returns the counts mapping to the bundle this transcript is in
     * the total bundle counts is the sum of this value for all transcripts in the bundle
     * @return a portion of the counts mapping to the bundle this transcript is in 
     */
    size_t bundle_counts() { return _bundle_counts; }
    
    /**
     * a member function that increases the expected fragment counts and variance by a given (logged) fragment mass
     * @param p a double for the (logged) probability that the fragment was generated by this transcript
     * @param mass a double specifying the (logged) mass of the fragment being mapped
     */
    void add_mass(double p, double mass);
    
    /**
     * a member function that increases the ccount of fragments uniquely mapped to this transcript
     * @param incr_amt a size_t to increase the counts by
     */
    void incr_uniq_counts(size_t incr_amt = 1)
    {
        _uniq_counts += incr_amt;
    }
    
    /**
     * a member function that increases the counts mapping to the bundle this transcript is in
     * the total bundle counts is the sum of this value for all transcripts in the bundle
     * @param incr_amt a size_t to increase the counts by
     */
    void incr_bundle_counts(size_t incr_amt = 1)
    {
        _bundle_counts += incr_amt;
    }
    
    /**
     * a member function that returns (a value proportional to) the log likelihood the given fragment
     * originated from this transcript
     * @param frag a FragHit to return the likelihood of being originated from this transcript
     * @return (a value proportional to) the log likelihood the given fragment originated from this transcript
     */
    double log_likelihood(const FragHit& frag) const;

    
    /**
     * a member function that calcualtes and returns the effective length of the transcript (non-logged)
     * @return the effective length of the transcript calculated as \f$ \tilde{l} = \sum_{l=1}^{L(t)}\sum_{i=1}^{L(t)} D(l)b_5[i]*b_3[i+l] \f$
     */
    double effective_length() const;
    
    /**
     * a member function that calcualtes and returns the estimated effective length of the transcript (non-logged) using the avg bias
     * @return the estimated effective length of the transcript calculated as \f$ \tilde{l} = \bar{bias}\sum_{l=1}^{L(t)} D(l)(L(t) - l + 1) \f$
     */
    double est_effective_length() const;
    
    /**
     * a member function that calcualtes and returns the effective length of the transcript (non-logged) ignoring bias and using the prior FLD distribution
     * @return the effective length of the transcript calculated as \f$ \tilde{l} = \sum_{l=1}^{L(t)} D_{prior}(l)(L(t) - l + 1) \f$
     */
    double unbiased_effective_length() const;

    /**
     * a member function that causes the transcript bias to be re-calculated by the _bias_table based on curent parameters
     */
    void update_transcript_bias();
    
};

typedef boost::unordered_map<TransID, Transcript*> TransMap;
typedef boost::unordered_map<size_t, double> TransPairMap;

// Typedefs to simplify the bundle partitioning
typedef boost::unordered_map<TransID, size_t> RankMap;
typedef boost::unordered_map<TransID, TransID> ParentMap;
typedef boost::associative_property_map<RankMap> Rank;
typedef boost::associative_property_map<ParentMap> Parent;
typedef boost::disjoint_sets<Rank, Parent> TransPartition;

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
     * a private map to look up pointers to Transcript objects by their hashed TransID id
     */
    TransMap _trans_map;
    
    /**
     * a private map to look up the covariance for pairs of Transcripts by their combined hashed TransIDs
     * these values are stored positive, even though they are negative
     */
    TransPairMap _covar_map;
    
    // objects used to partition the transcripts into bundles
    RankMap _rank_map;
    ParentMap _parent_map;
    Rank _rank;
    Parent _parent;
    TransPartition _bundles;
    
    /**
     * a private double specifying the initial Transcript mass pseudo-counts for each bp (non-logged)
     */
    double _alpha;
    
    /**
     * a private function that adds a transcript pointer to the table
     * @param trans a pointer to the transcript to be added
     */
    void add_trans(Transcript* trans);  
    
public:
    /**
     * TranscriptTable Constructor
     * @param trans_fasta_file a string storing the path to the fasta file from which to load transcripts
     * @param alpha a double that specifies the intial pseudo-counts for each bp of the transcripts (non-logged)
     * @param fld a pointer to the global Fragment Length Distribution (FLD) object
     * @param bias_table a pointer to the global BiasBoss object
     * @param mismatch_table a pointer to the global MismatchTable object
     */
    TranscriptTable(const std::string& trans_fasta_file, double alpha, const FLD* fld, BiasBoss* bias_table, const MismatchTable* mismatch_table);
    
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
     * a member function that returns the number of transcripts in the table
     * @return number of transcripts in the table
     */
    size_t size() const { return _trans_map.size(); }
    
    /**
     * a member function that increases the covariance between two transcripts by the specified (logged) amount
     * these values are stored positive even though they are negative
     * @param trans1 one of the transcripts in the pair
     * @param trans2 the other transcript in the pair
     * @param covar a double specifying the (logged) amount to increase the pair's covariance by
     */
    void update_covar(TransID trans1, TransID trans2, double covar);
    
    /**
     * a member function that returns the covariance between two transcripts
     * these returned value will be the log of the negative of the true value
     * @param trans1 one of the transcripts in the pair
     * @param trans2 the other transcript in the pair
     * @return the log of the negative of the pair's covariance
     */
    double get_covar(TransID trans1, TransID trans2);
    
    /**
     * a member function that returns the number of pairs of transcripts with non-zero covariance
     * @return the number of transcript pairs with non-zero covariance
     */
    size_t covar_size() const { return _covar_map.size(); }
    
    /**
     * a member function that returns the bundle representative of the given transcript in the partitioning
     * @param trans the TransID of the transcript whose representative is requested
     * @return the TransID of the representative for the bundle the given transcript is in
     */
    TransID get_trans_rep(TransID trans);
    
    /**
     * a member function that merges the bundles represented by the two given transcripts
     * @param rep1 the TransID of the first bundle representative
     * @param rep2 the TransID of the second bundle representative
     * @return the TransID of the representative for the new merged bundle
     */
    TransID merge_bundles(TransID rep1, TransID rep2);
    
    /**
     * a member function that returns the number of bundles in the partition
     * @return the number of bundles in the partition
     */
    size_t num_bundles();
    
    /**
     * a member function that outputs the bundles of the partition in a tab-delimited file called 'bundles.tab'
     * in the given output directory
     * each line contains a space-separated list of transcripts in a single bundle
     * @param output_dir the directory to output the bundle file to
     */
    void output_bundles(std::string output_dir);
    
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
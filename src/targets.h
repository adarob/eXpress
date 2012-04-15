//
//  targets.h
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
#include "sequence.h"

class FLD;
class FragHit;
class BiasBoss;
class MismatchTable;

/**
 *  The RoundParams struct stores the target parameters unique to a given round (iteration) of EM
 *  @author    Adam Roberts
 *  @date      2012
 *  @copyright Artistic License 2.0
 **/
struct RoundParams
{
    /**
     * a public double that stores the (logged) assigned mass based on observed fragment mapping probabilities
     */
    double mass;
        
    /**
     * a public double that stores the (logged) total mass of ambiguous fragments mapping to the target
     */
    double tot_ambig_mass;
   
    /**
     * a public double that stores the (logged) variance due to uncertainty on p
     */
    double mass_var;
    
    /**
     * a public double that stores the (logged) weighted sum of the variance on the assignments
     */
    double var_sum;
    
    /**
     * RoundParams constructor sets initial values for parameters
     */
    RoundParams() : mass(HUGE_VAL), tot_ambig_mass(HUGE_VAL), mass_var(HUGE_VAL), var_sum(HUGE_VAL) {}
};

typedef size_t TargID;


/**
 * The Target class is used to store objects for the targets
 * being mapped to.  Besides storing basic information about the object
 * (id, length), it also stores a mass based on the number of 
 * fragments mapping to the object.  To help with updating this number,
 * it returns the likelihood that a given fragment originated from it.
 * These values are stored and returned in log space.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class Target
{
    /**
     * a private pointer to the struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
     */
    const Globals* _globs;
    
    /**
     * a private TargID that stores the hashed target name
     */
    TargID _id;
    
    /**
     * a private string that stores the target name
     */
    std::string _name;
    
    /**
     * a private Sequence object that stores the forward target sequence
     */
    SequenceFwd _seq_f;

    /**
     * a private Sequence object that stores the reverse target sequence
     */
    SequenceRev _seq_r;
    
    /**
     * a private double object that stores the pseudo-mass-per-base
     */
    double _alpha;
    
    /**
     * a private RoundParams struct that stores the parameters for the current round
     */
    RoundParams _curr_params;
    
    /**
     * a private RoundParams struct that stores the parameters for the previous round
     */
    RoundParams _last_params;
    
    /**
     * a private pointer to the RoundParams that should be used in any accessors
     */
    RoundParams* _ret_params;
    
    /**
     * a private size_t that stores the number of fragments (non-logged) uniquely mapping
     * to this target
     */
    size_t _uniq_counts;
    
    /**
     * a private size_t that stores the fragment counts (non-logged) for the bundle
     * the total bundle counts is the sum of this value for all targets in the bundle
     */
    size_t _tot_counts;
    
    /**
     * a private pointer to the Bundle this Target is a member of
     */
    Bundle* _bundle;
    
    /**
     * a private mutex to provide thread-safety for variables with threaded update
     */
    mutable boost::mutex _lock;
    
    /**
     * a point to a private double vector storing the (logged) 5' bias at each position (accessed by multiple threads)
     */
    std::vector<float>* _start_bias;
    
    /**
     * a pointer to a private double vector storing the (logged) 3' bias at each position (accessed by multiple threads)
     */
    std::vector<float>* _end_bias;
    
    /**
     * a private double storing the (logged) product of the average 3' and 5' biases for the target (accessed by multiple threads)
     */
    double _avg_bias;
    
    /**
     * a private double storing the most recently updated (logged) effective length as calculated by the bias updater thread
     */
    double _cached_eff_len;
    
    /**
     * a private boolean specifying whether a unique solution exists, always false unless a unique read is mapped to the target
     * or all other targets in a mapping are solveable
     */
    bool _solveable;
    
public:
    
    /**
     * Target Constructor
     * @param name a string that stores the target name
     * @param seq a string that stores the target sequence
     * @param DOC
     * @param alpha a double that specifies the intial pseudo-counts (non-logged)
     * @param globs a pointer to the struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
     */
    Target(const size_t id, const std::string& name, const std::string& seq, bool prob_seq, double alpha, const Globals* globs);
    
    /**
     * Target Destructor deletes bias vectors
     */    
    ~Target()
    {
        if (_start_bias)
            delete _start_bias;
        if (_end_bias)
            delete _end_bias;
    }
    
    /**
     * a member function that locks the target mutex to provide thread safety
     */
    void lock() { _lock.lock(); }
    
    /**
     * a member function that unlocks the target mutex to provide thread safety
     */
    void unlock() { _lock.unlock(); }

    /**
     * a member function that returns the target name
     * @return string containing target name
     */
    const std::string& name() const { return _name; }
    
    /**
     * a member function that returns the target id
     * @return TargID target ID
     */
    TargID id() const { return _id; }
    
    /**
     * a member function that returns the target sequence (const)
     * @return string containing target sequence
     */
    const Sequence& seq(bool rev=false) const { if (rev) return _seq_r; return _seq_f; }
    
    /**
     * a member function that returns the target sequence (non-const)
     * @return string containing target sequence
     */
    Sequence& seq(bool rev) { if (rev) return _seq_r; return _seq_f; }
    
    /**
     * a member function that returns the target length
     * @return target length
     */
    size_t length() const { return _seq_f.length(); }
    
    /**
     * a member function that returns the current estimated rho (logged, w/ pseudo-counts) for the target
     * @return the current estimated rho
     */
    double rho() const;
    
    /**
     * a member function that returns the current (logged) probabilistically assigned fragment mass
     * @param with_pseudo a boolean specifying whether pseudo-counts should be included in returned mass
     * @return logged mass
     */
    double mass(bool with_pseudo=true) const;
    
    /**
     * a member function that returns the total (logged) variance on mass
     * @return the total (logged) variance on mass
     */
    double mass_var(bool with_pseudo=true) const;
    
    /**
     * a member function that returns the (logged) weighted sum of the variance on the assignments
     * @return the (logged) weighted sum of the variance on the assignments
     */
    double var_sum() const { return _ret_params->var_sum; }
            
    /**
     * a member function that returns the (logged) total mass of ambiguous fragments mapping to the target
     * @return the (logged) total mass of ambiguous fragments mapping to the target
     */
    double tot_ambig_mass() const { return _ret_params->tot_ambig_mass; }
    
    /**
     * a member function that readies the target object for the next round of batch EM
     */
    void round_reset();
    
    /**
     * a member function that returns the current count of fragments mapped to this target (uniquely or ambiguously)
     * @return total fragment count
     */
    size_t tot_counts() const { return _tot_counts; }
    
    /**
     * a member function that returns the current count of fragments uniquely mapped to this target
     * @return unique fragment count
     */
    size_t uniq_counts() const { return _uniq_counts; }
    
    /**
     * a member function that returns the Bundle this Target is a member of
     * @return a pointer to the Bundle this target is a member of
     */
    Bundle* bundle() { return _bundle; }

    /**
     * a member function that set the Bundle this Target is a member of
     * @param b a pointer to the Bundle to set this Target as a member of
     */
    void bundle(Bundle* b) { _bundle = b; }
    
    /**
     * a member function that increases the expected fragment counts and variance by a given (logged) fragment mass
     * @param p a double for the (logged) probability that the fragment was generated by this target
     * @param v a double for the (logged) approximate variance (uncertainty) on the probability
     * @param mass a double specifying the (logged) mass of the fragment being mapped
     */
    void add_mass(double p, double v, double mass);
    
    /**
     * a member function that increases the count of fragments mapped to this target
     * @param uniq a bool specifying whether or not the fragment uniquely maps to this target
     * @param incr_amt a size_t to increase the counts by
     */
    void incr_counts(bool uniq, size_t incr_amt = 1)
    {
        if (uniq)
            _solveable = true;
        _tot_counts += incr_amt;
        _uniq_counts += incr_amt*uniq;
    }
    
    /**
     * a member function that returns (a value proportional to) the log likelihood the given fragment
     * originated from this target
     * @param frag a FragHit to return the likelihood of being originated from this target
     * @param with_pseudo a FragHit specifying whether or not pseudo-counts should be included in the calculation
     * @return (a value proportional to) the log likelihood the given fragment originated from this target
     */
    double log_likelihood(const FragHit& frag, bool with_pseudo) const;
        
    /**
     * a member function that calcualtes and returns the estimated effective length of the target (logged) using the avg bias
     * @param fld an optional pointer to a different FLD than the global one, for thread-safety
     * @param with_bias a boolean specifying whether or not the average bias should be used in the calculation
     * @return the estimated effective length of the target calculated as \f$ \tilde{l} = \bar{bias}\sum_{l=1}^{L(t)} D(l)(L(t) - l + 1) \f$
     */
    double est_effective_length(FLD* fld = NULL, bool with_bias=true) const;
    
    /**
     * a member function that returns the most recently estimated effective length (logged) as calculated by the bias updater thread 
     * @return the cached effective length of the target calculated
     */
    double cached_effective_length(bool with_bias=true) const;

    /**
     * a member function that causes the target bias to be re-calculated by the _bias_table based on curent parameters
     * @param bias_table an optional pointer to a different BiasBoss than the global one, for thread-safety
     * @param fld an optional pointer to a different FLD than the global one, for thread-safety
     */
    void update_target_bias(BiasBoss* bias_table = NULL, FLD* fld = NULL);

    /**
     * a member function that returns the _solveable flag
     * @return a boolean specifying whether or not the target has a unique solution
     */    
    bool solveable() { return _solveable; }
    
    /**
     * a member function that sets the _solveable flag
     * @param a boolean specifying whether or not the target has a unique solution
     */   
    void solveable(bool s) { _solveable = s; }
    
};

typedef std::vector<Target*> TransMap;
typedef boost::unordered_map<std::string, size_t> TransIndex;
typedef boost::unordered_map<size_t, float> CovarMap;
typedef boost::unordered_map<std::string, double> AlphaMap;

/**
 * The TargetTable class is used to keep track of the Target objects for a run.
 * The constructor parses a fasta file to generate the Target objects and store them in a map
 * that allows them to be looked up based on their string id.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class TargetTable
{
    /**
     * a private pointer to the struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
     */
    Globals* _globs;
        
    /**
     * a private map to look up pointers to Target objects by their TargID id
     */
    TransMap _targ_map;
    
    /**
     * the private table to keep track of Bundles
     */
    BundleTable _bundle_table;
    
    /**
     * a private CovarTable to look up the covariance for pairs of Targets by their combined hashed TargIDs
     * these values are stored positive, even though they are negative
     */
    CovarTable _covar_table;
    
    /**
     * a private double that stores the (logged) total mass per base (including pseudo-counts) to allow for rho calculations
     */
    double _total_fpb;

    /**
     * a private mutex to make accesses to _total_fpb thread-safe
     */
    mutable boost::mutex _fpb_mut;
    
    /**
     * a private function that validates and adds a target pointer to the table
     * @param name the name of the trancript
     * @param seq the sequence of the target
     * @param DOC
     * @param alpha a double that specifies the initial pseudo-counts for each bp of the target (non-logged)
     * @param targ_index the target-to-index map from the alignment file
     * @param targ_lengths the target-to-length map from the alignment file (for validation)
     */
    void add_targ(const std::string& name, const std::string& seq, bool prob_seqs, double alpha, const TransIndex& targ_index, const TransIndex& targ_lengths);  
    
public:
    /**
     * TargetTable Constructor
     * @param targ_fasta_file a string storing the path to the fasta file from which to load targets
     * @param targ_index the target-to-index map from the alignment file
     * @param targ_lengths the target-to-length map from the alignment file
     * @param DOC
     * @param alpha a double that specifies the intial pseudo-counts for each bp of the targets (non-logged)
     * @param alpha_map an optional pointer to a map object that specifies proportional weights of pseudo-counts for each target
     * @param globs a pointer to the struct containing pointers to the global parameter tables (bias_table, mismatch_table, fld)
     */
    TargetTable(const std::string& targ_fasta_file, const TransIndex& targ_index, const TransIndex& targ_lengths, bool prob_seqs, double alpha, const AlphaMap* alpha_map, Globals* globs);
    
    /**
     * TargetTable Destructor
     * deletes all of the target objects in the table
     */
    ~TargetTable();
        
    /**
     * a member function that returns a pointer to the target with the given id 
     * @param id of the target queried
     * @return pointer to the target wit the given id
     */
    Target* get_targ(TargID id);
    
    /**
     * a member function that readies all Target objects in the table for the next round of batch EM
     */
    void round_reset();
    
    /**
     * a member function that returns the number of targets in the table
     * @return number of targets in the table
     */
    size_t size() const { return _targ_map.size(); }
    
    /**
     * a member function that returns the (logged) total mass per base (including pseudo-counts)
     * @return the (logged) total mass per base (including pseudo-counts)
     */    
    double total_fpb() const;
    
    /**
     * a member function that increments the (logged) total mass per base
     * @param incr_amt the (logged) amount to increment by
     */   
    void update_total_fpb(double incr_amt);
    
    /**
     * a member function that increases the covariance between two targets by the specified amount
     * these values are stored positive even though they are negative (logged)
     * @param targ1 one of the targets in the pair
     * @param targ2 the other target in the pair
     * @param covar a double specifying the amount to increase the pair's covariance by (logged)
     */
    void update_covar(TargID targ1, TargID targ2, double covar) { _covar_table.increment(targ1, targ2, covar); }
    
    /**
     * a member function that returns the covariance between two targets
     * these returned value will be the log of the negative of the true value
     * @param targ1 one of the targets in the pair
     * @param targ2 the other target in the pair
     * @return the negative of the pair's covariance (logged)
     */
    double get_covar(TargID targ1, TargID targ2) { return _covar_table.get(targ1, targ2); }
    
    /**
     * a member function that returns the number of pairs of targets with non-zero covariance
     * @return the number of target pairs with non-zero covariance
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
     * a member function that outputs the final expression data in a file called 'results.xprs'
     * and (optionally) the variance-covariance matrix in 'varcov.xprs' in the given output directory
     * @param output_dir the directory to output the expression file to
     * @param tot_counts the total number of observed mapped fragments
     * @param output_varcov boolean specifying whether to also output the variance-covariance matrix
     */
    void output_results(std::string output_dir, size_t tot_counts, bool output_varcov);
    
    /**
     * a member function for driving a thread that continuously updates the target bias values
     */
    void threaded_bias_update(boost::mutex* mut);
};

#endif
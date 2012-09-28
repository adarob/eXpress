/**
 *  targets.h
 *  express
 *
 *  Created by Adam Roberts on 3/20/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 */

#ifndef TRANSCRIPTS_H
#define TRANSCRIPTS_H

#include <boost/scoped_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/unordered_map.hpp>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "main.h"
#include "bundles.h"
#include "sequence.h"

class FLD;
struct FragHit;
class BiasBoss;
class MismatchTable;
class Librarian;

/**
 * The RoundParams struct stores the target parameters unique to a given round
 * (iteration) of EM.
 * @author    Adam Roberts
 * @date      2012
 * @copyright Artistic License 2.0
 **/
struct RoundParams {
  /**
   * A public double that stores the (logged) assigned mass based on observed
   * fragment mapping probabilities.
   */
  double mass;
  /**
   * A public double that stores the (logged) total mass of ambiguous fragments
   * mapping to the target.
   */
  double tot_ambig_mass;
  /**
   * A public double that stores the (logged) variance due to uncertainty on p.
   */
  double mass_var;
  /**
   * A public double that stores the (logged) weighted sum of the variance on
   * the assignments.
   */
  double var_sum;
  /**
   * RoundParams constructor sets initial values for parameters
   */
  RoundParams() : mass(LOG_0), tot_ambig_mass(LOG_0), mass_var(LOG_0),
                  var_sum(LOG_0) {}
};

typedef size_t TargID;

/**
 * The Target class is used to store objects for the targets being mapped to.
 * Besides storing basic information about the object (id, length), it also
 * stores a mass based on the number of fragments mapping to the object as well
 * as parameters for variance. To help with updating these values, it computes
 * the likelihood that a given fragment originated from it. These values are
 * stored and returned in log space.
 *  @author  Adam Roberts
 *  @date    2011
 *  @copyright Artistic License 2.0
 **/
class Target {
  /**
   * A private pointer to the struct containing pointers to the global
   * parameter tables (bias_table, mismatch_table, fld).
   */
  const Librarian* _libs;
  /**
   * A private TargID that stores the hashed target name.
   */
  TargID _id;
  /**
   * A private string that stores the target name.
   */
  std::string _name;
  /**
   * A private Sequence object that stores the forward target sequence.
   */
  SequenceFwd _seq_f;
  /**
   * A private Sequence object that stores the reverse target sequence.
   */
  SequenceRev _seq_r;
  /**
   * A private double object that stores the pseudo-mass-per-base.
   */
  double _alpha;
  /**
   * A private RoundParams struct that stores the parameters for the current
   * round.
   */
  RoundParams _curr_params;
  /**
   * A private RoundParams struct that stores the parameters for the previous
   * round.
   */
  RoundParams _last_params;
  /**
   * A private pointer to the RoundParams that should be used in any
   * accessors.
   */
  RoundParams* _ret_params;
  /**
   * A private size_t that stores the number of fragments (non-logged)
   * uniquely mapping to this target.
   */
  size_t _uniq_counts;
  /**
   * A private size_t that stores the fragment counts (non-logged) for the
   * bundle. The total bundle counts is the sum of this value for all targets
   * in the bundle.
   */
  size_t _tot_counts;
  /**
   * A private pointer to the Bundle this Target is a member of.
   */
  Bundle* _bundle;
  /**
   * A private mutex to provide thread-safety for variables with threaded
   * update.
   */
  mutable boost::mutex _mutex;
  /**
   * A scoped pointer to a private float vector storing the (logged) 5' bias
   * at each position.
   */
  boost::scoped_ptr<std::vector<float> > _start_bias;
  /**
   * A scoped pointer to a private float vector storing the (logged) 3' bias
   * at each position.
   */
  boost::scoped_ptr<std::vector<float> > _end_bias;
  /**
   * A private double storing the (logged) product of the average 3' and 5'
   * biases for the target.
   */
  double _avg_bias;
  /**
   * A private double storing the most recently updated (logged) effective
   * length as calculated by the bias updater thread.
   */
  double _cached_eff_len;
  /**
   * A private boolean specifying whether a unique solution exists. True iff
   * a unique read is mapped to the target or all other targets in a mapping
   * are solvable.
   */
  bool _solvable;

public:
  /**
   * Target Constructor.
   * @param id a unique TargID identifier.
   * @param name a string that stores the target name.
   * @param seq a string that stores the target sequence.
   * @param prob_seq a bool that specifies if the sequence is to be treated
   *        probablistically. For RDD detection.
   * @param alpha a double that specifies the intial pseudo-counts
   *        (non-logged).
   * @param libs a pointer to the struct containing pointers to the global
   *        parameter tables (bias_table, mismatch_table, fld).
   */
  Target(TargID id, const std::string& name, const std::string& seq,
       bool prob_seq, double alpha, const Librarian* libs);
  /**
   * A member function that locks the target mutex to provide thread safety.
   * The lock should be held by any thread that calls a method of the Target.
   */
  void lock() { _mutex.lock(); }
  /**
   * A member function that unlocks the target mutex.
   */
  void unlock() { _mutex.unlock(); }
  /**
   * An accessor for the target name.
   * @return string containing target name.
   */
  const std::string& name() const { return _name; }
  /**
   * An accessor for the target id.
   * @return The target ID.
   */
  TargID id() const { return _id; }
  /**
   * An accessor for the the target's Sequence (const).
   * @param rev a bool specifying whether to return the reverse complement.
   * @return Const reference to the target's Sequence object.
   */
  const Sequence& seq(bool rev=false) const {
    if (rev) {
      return _seq_r;
    }
    return _seq_f;
  }
  /**
   * An accessor for the the target's Sequence (non-const).
   * @param rev a bool specifying whether to return the reverse complement.
   * @return Non-const reference to the target's Sequence object.
   */
  Sequence& seq(bool rev) {
    if (rev) {
      return _seq_r;
    }
    return _seq_f;
  }
  /**
   * An accessor for the length of the target sequence.
   * @return The target sequence length.
   */
  size_t length() const { return _seq_f.length(); }
  /**
   * An accessor for he current estimated relative abundance of the target. The
   * value is logged and includes pseudo-counts.
   * @return The current estimated rho.
   */
  double rho() const;
  /**
   * An accessor for the current (logged) probabilistically assigned fragment
   * mass.
   * @param with_pseudo a boolean specifying whether pseudo-counts should be
   *        included in returned mass.
   * @return The logged mass.
   */
  double mass(bool with_pseudo=true) const;
  /**
   * An accessor for the total (logged) variance on mass.
   * @param with_pseudo a boolean specifying whether pseudo-counts should be
   *        included in returned variance.
   * @return The total (logged) variance on mass.
   */
  double mass_var(bool with_pseudo=true) const;
  /**
   * An accessor for the (logged) weighted sum of the variance on assignments.
   * @return The (logged) weighted sum of the variance on the assignments.
   */
  double var_sum() const { return _ret_params->var_sum; }
  /**
   * An accessor for the the (logged) total mass of ambiguous fragments mapping
   * to the target.
   * @return The (logged) total mass of ambiguous fragments mapping to the
   *         target.
   */
  double tot_ambig_mass() const { return _ret_params->tot_ambig_mass; }
  /**
   * A member function that prepares the target object for the next round of
   * batch EM.
   */
  void round_reset();
  /**
   * An accessor for the current count of fragments mapped to this target
   * either uniquely or ambiguously.
   * @return The total fragment count.
   */
  size_t tot_counts() const { return _tot_counts; }
  /**
   * An accessor for the the current count of fragments uniquely mapped to this
   * target.
   * @return The unique fragment count.
   */
  size_t uniq_counts() const { return _uniq_counts; }
  /**
   * An accessor for the pointer to the Bundle this Target is a member of.
   * @return A pointer to the Bundle this target is a member of.
   */
  Bundle* bundle() const { return _bundle; }
  /**
   * A mutator to set the Bundle this Target is a member of.
   * @param b a pointer to the Bundle to set this Target as a member of.
   */
  void bundle(Bundle* b) { _bundle = b; }
  /**
   * A member function that increases the expected fragment counts and
   * variance by a given (logged) fragment mass.
   * @param p a double for the (logged) probability that the fragment was
   *        generated by this target.
   * @param v a double for the (logged) approximate variance (uncertainty) on
   *        the probability p.
   * @param mass a double specifying the (logged) mass of the fragment being
   *        mapped.
   */
  void add_mass(double p, double v, double mass);
  /**
   * A member function that increases the count of fragments mapped to this
   * target.
   * @param uniq a bool specifying whether or not the fragment uniquely maps
   *        to this target.
   * @param incr_amt a size_t to increase the counts by.
   */
  void incr_counts(bool uniq, size_t incr_amt = 1) {
    if (uniq) {
      _solvable = true;
    }
    _tot_counts += incr_amt;
    _uniq_counts += incr_amt * uniq;
  }
  /**
   * A member function that returns (a value proportional to) the log likelihood
   * the given fragment originated from this target.
   * @param frag a FragHit to return the likelihood of originating from this
   *        target.
   * @param with_pseudo a bool specifying whether or not pseudo-counts should be
   *        included in the calculation.
   * @return A value proportional to the log likelihood the given fragment
   *         originated from this target.
   */
  double log_likelihood(const FragHit& frag, bool with_pseudo) const;
  /**
   * A member function that calculates and returns the estimated effective
   * length of the target (logged) using the average bias.
   * @param fld an optional pointer to a different FLD than the global one, for
   *        thread-safety.
   * @param with_bias a boolean specifying whether or not the average bias
   *        should be included in the return value.
   * @return The estimated effective length of the target calculated as
   *         \f$ \tilde{l} = \bar{bias}\sum_{l=1}^{L(t)} D(l)(L(t) - l + 1) \f$.
   */
  double est_effective_length(FLD* fld = NULL, bool with_bias=true) const;
  /**
   * An accessor for the most recently estimated effective length (logged) as
   * calculated by the bias updater thread.
   * @param with_bias a boolean specifying whether or not the average bias
   * should be included in the return value.
   * @return The cached effective length of the target.
   */
  double cached_effective_length(bool with_bias=true) const;
  /**
   * A member function that causes the target bias to be re-calculated by the
   * _bias_table based on curent parameters.
   * @param bias_table an optional pointer to a different BiasBoss than the
   *        global one, for thread-safety.
   * @param fld an optional pointer to a different FLD than the global one,
   *        for thread-safety.
   */
  void update_target_bias(BiasBoss* bias_table = NULL, FLD* fld = NULL);
  /**
   * An accessor for the _solvable flag.
   * @return a boolean specifying whether or not the target has a unique
   *         solution for its abundance estimate.
   */

  bool solvable() { return _solvable; }
  /**
   * A mutator that sets the _solvable flag.
   * @param a boolean specifying whether or not the target has a unique solution
   *        for its abundance estimate.
   */

  void solvable(bool s) { _solvable = s; }
};

typedef std::vector<Target*> TransMap;
typedef boost::unordered_map<std::string, size_t> TransIndex;
typedef boost::unordered_map<size_t, float> CovarMap;
typedef boost::unordered_map<std::string, double> AlphaMap;

/**
 * The TargetTable class is used to keep track of the Target objects for a run.
 * The constructor parses a fasta file to generate the Target objects and stores
 * them in a map keyed by their string id.
 *  @author  Adam Roberts
 *  @date    2011
 *  @copyright Artistic License 2.0
 **/
class TargetTable {
  /**
   * A private pointer to the struct containing pointers to the global parameter
   * tables (bias_table, mismatch_table, fld).
   */
  const Librarian* _libs;
  /**
   * A private map to look up pointers to Target objects by their TargID id.
   */
  TransMap _targ_map;
  /**
   * The private table to keep track of Bundle objects.
   */
  BundleTable _bundle_table;
  /**
   * A private table to look up the covariance for pairs of Targets by their
   * combined hashed TargIDs. These values are stored logged and positive, even
   * though they are negative.
   */
  CovarTable _covar_table;
  /**
   * A private double that stores the (logged) total mass per base
   * (including pseudo-counts) to allow for rho calculations.
   */
  double _total_fpb;
  /**
   * A private mutex to make accesses to _total_fpb thread-safe.
   */
  mutable boost::mutex _fpb_mut;

  /**
   * A private function that validates and adds a target pointer to the table.
   * @param name the name of the trancript.
   * @param seq the sequence of the target.
   * @param prob_seqs a bool that specifies if the sequence is to be treated
   *        probablistically, for RDD detection.
   * @param alpha a double that specifies the initial pseudo-counts for each bp
   *        of the target (non-logged).
   * @param targ_index the target-to-index map from the alignment file.
   * @param targ_lengths the target-to-length map from the alignment file, for
   *        validation.
   */
  void add_targ(const std::string& name, const std::string& seq, bool prob_seqs,
                double alpha, const TransIndex& targ_index,
                const TransIndex& targ_lengths);

public:
  /**
   * TargetTable Constructor.
   * @param targ_fasta_file a string storing the path to the fasta file from
   *        which to load targets.
   * @param prob_seqs a bool that specifies if the sequence is to be treated
   *        probablistically, for RDD detection.
   * @param alpha a double that specifies the intial pseudo-counts for each bp
   *        of the targets (non-logged).
   * @param alpha_map an optional pointer to a map object that specifies
   *        proportional weights of pseudo-counts for each target.
   * @param libs a pointer to the struct containing pointers to the global
   *        parameter tables (bias_table, mismatch_table, fld).
   */
  TargetTable(const std::string& targ_fasta_file, bool prob_seqs, double alpha,
              const AlphaMap* alpha_map, const Librarian* libs);
  /**
   * TargetTable Destructor. Deletes all of the target objects in the table.
   */
  ~TargetTable();
  /**
   * A member function that returns a pointer to the target with the given id.
   * @param id of the target queried.
   * @return A pointer to the target with the given id.
   */
  Target* get_targ(TargID id);
  /**
   * A member function that readies all Target objects in the table for the next
   * round of batch EM.
   */
  void round_reset();
  /**
   * An accessor for the number of targets in the table.
   * @return The number of targets in the table.
   */
  size_t size() const { return _targ_map.size(); }
  /**
   * An accessor for the (logged) total mass per base, including pseudo-counts.
   * @return The (logged) total mass per base, including pseudo-counts.
   */

  double total_fpb() const;
  /**
   * a member function that increments the (logged) total mass per base.
   * @param incr_amt the (logged) amount to increment by.
   */

  void update_total_fpb(double incr_amt);
  /**
   * A member function that increases the (logged) covariance between two
   * targets by the specified amount. These values are stored positive even
   * though they are negative.
   * @param targ1 one of the targets in the pair
   * @param targ2 the other target in the pair
   * @param covar a double specifying the amount to increase the pair's
   *        covariance by (logged)
   */
  void update_covar(TargID targ1, TargID targ2, double covar) {
    _covar_table.increment(targ1, targ2, covar);
  }
  /**
   * An accessor for the covariance between two targets. These returned value
   * will be the log of the negative of the true value.
   * @param targ1 one of the targets in the pair.
   * @param targ2 the other target in the pair.
   * @return The negative of the pair's covariance (logged).
   */
  double get_covar(TargID targ1, TargID targ2) {
    return _covar_table.get(targ1, targ2);
  }
  /**
   * An accessor for number of pairs of targets with non-zero covariance.
   * @return The number of target pairs with non-zero covariance.
   */
  size_t covar_size() const { return _covar_table.size(); }
  /**
   * A member function that merges the given Bundles.
   * @param b1 a pointer to the first Bundle to merge.
   * @param b2 a pointer to the second Bundle to merge.
   * @return A pointer to the merged Bundle.
   */
  Bundle* merge_bundles(Bundle* b1, Bundle* b2);
  /**
   * An accessor for the number of bundles in the partition.
   * @return The number of bundles in the partition.
   */
  size_t num_bundles() const { return _bundle_table.size(); }
  /**
   * A member function that outputs the final expression data in a file called
   * 'results.xprs', (optionally) the variance-covariance matrix in
   * 'varcov.xprs', and (optionally) the RDD p-values in the given output
   * directory.
   * @param output_dir the directory to output the expression file to.
   * @param tot_counts the total number of observed mapped fragments.
   * @param output_varcov boolean specifying whether to also output the
   *        variance-covariance matrix
   * @param output_rdds boolean specifying whether to also output the RDD
   *        p-values.
   */
  void output_results(std::string output_dir, size_t tot_counts,
                      bool output_varcov=false, bool output_rdds=false);
  /**
   * A member function to be run asynchronously that continuously updates the
   * background bias values, target bias values, and target effective lengths.
   * @param mutex a pointer to the mutex to be used to protect the global fld
   *        and bias tables during updates.
   */
  void asynch_bias_update(boost::mutex* mutex);
};

#endif

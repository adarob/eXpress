/**
 *  sequence.h
 *  express
 *
 *  Created by Adam Roberts on 1/24/12.
 *  Copyright 2012 Adam Roberts. All rights reserved.
 **/

#ifndef express_sequence_h
#define express_sequence_h

#include <boost/scoped_array.hpp>
#include <string>
#include "frequencymatrix.h"

/**
 * Helper function to encode a nucleotide character to a size_t value.
 * @param c the nucleotide character to be encoded.
 * @return A char value encoding the nucleotide.
 */
inline char ctoi(const char c) {
  switch(c) {
    case 'A':
    case 'a':
       return 0;
    case 'C':
    case 'c':
       return 1;
    case 'G':
    case 'g':
       return 2;
    case 'T':
    case 't':
       return 3;
    default:
       return 0;
  }
}
/**
 * Helper function to return the encoded complement of the given encoded
 * nucleotide.
 * @param c the encoded nucleotide to complement.
 * @return The encoded complement of the given encoded nucleotide.
 */
inline char complement(const char c) {
  return c^3;
}

/**
 * The Sequence class is an abstract class whose implmentations are used to
 * store and access encoded nucleotide sequences. They also supports
 * probabilistic sequences, meaning that each position stores a distribution
 * over nucleotides.
 *  @author    Adam Roberts
 *  @date      2012
 *  @copyright Artistic License 2.0
 **/
class Sequence {
 public:
  /**
   * Dummy Sequence constructor.
   */
  Sequence() {}
  /**
   * Dummy Sequence destructor.
   */
  virtual ~Sequence() {}
  /**
   * A member function that encodes the given sequence and overwrites the
   * current stored sequence with it.
   * @param seq the nucleotide sequence to encode and store.
   * @param rev a boolean if the sequence should be reverse complemented before
   *        encoding.
   */
  virtual void set(const std::string& seq, bool rev) = 0;
  /**
   * An accessor for the encoded character at the given index. If the sequence
   * is probabilistic, returns the mode. Otherwise, returns the reference.
   * @param index the index of the encoded character to return (assumed to be
   *        < _len)
   * @return The encoded character at the given index.
   */
  virtual size_t operator[](const size_t index) const = 0;
  /**
   * An accessor for the encoded reference character at the given index. May
   * differ from the operator[] if the sequence is probabilistic.
   * @param index the index of the encoded reference character to return
   *        (assumed to be < _len)
   * @return The encoded reference character at the given index.
   */
  virtual size_t get_ref(const size_t index) const = 0;
  /**
   * A member function that updates the posterior nucleotide distribution if
   * probabilistic.
   * @param index the index of the position to update.
   * @param nuc the nucleotide to increment the mass of.
   * @param mass the amount to increment the nucleotide's mass by.
   */
  virtual void update_est(const size_t index, const size_t nuc, float mass) = 0;
  /**
   * A member function that updates the observed frequency distribution if
   * probabilistic.
   * @param index the index of the position to update.
   * @param nuc the nucleotide to increment the mass of.
   * @param mass the amount to increment the nucleotide's mass by.
   */
  virtual void update_obs(const size_t index, const size_t nuc, float mass) = 0;
  /**
   * A member function that updates the expected frequency distribution if
   * probabilistic.
   * @param index the index of the position to update.
   * @param nuc the nucleotide to increment the mass of.
   * @param mass the amount to increment the nucleotide's mass by.
   */
  virtual void update_exp(const size_t index, const size_t nuc, float mass) = 0;
  /**
   * An accessor for the posterior nucleotide distribution (logged).
   * @param index the index of the position to access.
   * @param nuc the nucleotide to return the probability of.
   * @return The logged posterior probability of the given nucleotide at the
   *         given position.
   */
  virtual float get_prob(const size_t index, const size_t nuc) const = 0;
  /**
   * An accessor for the observed nucleotide frequency (logged).
   * @param index the index of the position to access.
   * @param nuc the nucleotide to return the frequency of.
   * @return The logged observed frequency of the given nucleotide a the given
   *         position.
   */
  virtual float get_obs(const size_t index, const size_t nuc) const = 0;
  /**
   * An accessor for the expected nucleotide frequency (logged).
   * @param index the index of the position to access.
   * @param nuc the nucleotide to return the frequency of.
   * @return The logged expected frequency of the given nucleotide a the given
   *         position.
   */
  virtual float get_exp(const size_t index, const size_t nuc) const = 0;
  /**
   * Accessor to determine if the sequence is probabilistic.
   * @return True iff the sequence is probabilistic.
   */
  virtual bool prob() const = 0;
  /**
   * An accessor for the length of the encoded sequence.
   * @return The length of the encoded sequence.
   */
  virtual size_t length() const = 0;
  /**
   * Accessor to determine if the sequence has 0 length.
   * @return True iff the sequence has 0 length.
   */
  virtual bool empty() const = 0;
  /**
   * A member function that calculates p-values based on the observed and
   * expected nucleotide frequences for the sequence. Experimental.
   * @param p_vals a reference to an empty vector of doubles to fill with
   *        p-values at each position.
   */
  virtual void calc_p_vals(std::vector<double>& p_vals) const = 0;
};

/**
 * The SequenceFwd class implements the Sequence abstract class for storing the
 * forward sequence. Documentation is only provided for methods not documented
 * in the abstract Sequence class.
 *  @author    Adam Roberts
 *  @date      2012
 *  @copyright Artistic License 2.0
 **/
class SequenceFwd: public Sequence
{
  /**
   * A char array that stores the encoded sequence (not null terminated).
   * Deleted with this.
   */
  boost::scoped_array<const char> _ref_seq;
  /**
   * A private bool specifying if the sequence is probabilistic (true) or fixed
   * to the reference (false).
   */
  bool _prob;
  /**
   * A private FrequencyMatrix to store the posterior nucleotide distributions
   * (if _prob).
   */
  FrequencyMatrix<float> _est_seq;
  /**
   * A private FrequencyMatrix to store the observed nucleotide frequencies
   * (if _prob).
   */
  FrequencyMatrix<float> _obs_seq;
  /**
   * A private FrequencyMatrix to store the expected nucleotide frequencies
   * (if _prob).
   */
  FrequencyMatrix<float> _exp_seq;
  /**
   * A private size_t storing the number of nucleotides in the sequence.
   */
  size_t _len;
    
 public:
  /**
   * Dummy SequenceFwd constructor.
   */
  SequenceFwd();
  /**
   * SequenceFwd constructor encodes and stores the given nucleotide sequence.
   * @param seq the nucleotide sequence string to encode and store.
   * @param rev a boolean if the sequence should be reverse complemented
   *        before encoding.
   */
  SequenceFwd(const std::string& seq, bool rev, bool prob=false);
  /**
   * SequenceFwd copy constructor.
   * @param other the SequenceFwd object to copy.
   */
  SequenceFwd(const SequenceFwd& other);
  /**
   * SequenceFwd assignment constructor, copies the given SequenceFwd object.
   * @param other the Sequence object to copy.
   */
  SequenceFwd& operator=(const SequenceFwd& other);
  // The following methods are documented in the abstract Sequence class.
  void set(const std::string& seq, bool rev);
  size_t operator[](const size_t index) const;
  size_t get_ref(const size_t index) const;
  float get_exp(const size_t index, const size_t nuc) const;
  float get_obs(const size_t index, const size_t nuc) const;
  void update_est(const size_t index, const size_t nuc, float mass);
  void update_obs(const size_t index, const size_t nuc, float mass);
  void update_exp(const size_t index, const size_t nuc, float mass);
  float get_prob(const size_t index, const size_t nuc) const;
  bool prob() const { return _prob; }
  bool empty() const { return _len==0; }
  size_t length() const { return _len; }
  void calc_p_vals(std::vector<double>& p_vals) const;
};

/**
 * The SequenceRev class implements the Sequence abstract class for accessing
 * the reverse sequence. Documentation is only provided for methods not
 * documented in the abstract Sequence class and SequenceFwd class. This class
 * acts by storing a pointer to a SequenceFwd object and reverse complementing 
 * the input and output appropriately.
 *  @author    Adam Roberts
 *  @date      2012
 *  @copyright Artistic License 2.0
 **/
class SequenceRev: public Sequence {
  /**
   * A private pointer to the SequenceFwd object this object is the reverse
   * complement of.
   */
  SequenceFwd* _seq;
    
 public:
  SequenceRev() {}
  SequenceRev(SequenceFwd& seq) : _seq(&seq) {}
  size_t length() const { return _seq->length(); }
  bool empty() const { return _seq->empty(); }
  // Set is not allowed in this class.
  void set(const std::string& seq, bool rev) { assert(false); exit(1); }
  size_t operator[](const size_t index) const {
    return complement(_seq->operator[](length()-index-1)); }
  size_t get_ref(const size_t index) const {
    return complement(_seq->get_ref(length()-index-1)); }
  float get_obs(const size_t index, const size_t nuc) const {
    return _seq->get_obs(length()-index-1, complement(nuc)); }
  float get_exp(const size_t index, const size_t nuc) const {
    return _seq->get_exp(length()-index-1, complement(nuc)); }
  void update_est(const size_t index, const size_t nuc, float mass) {
    _seq->update_est(length()-index-1, complement(nuc), mass); }
  void update_obs(const size_t index, const size_t nuc, float mass) {
    _seq->update_obs(length()-index-1, complement(nuc), mass); }
  void update_exp(const size_t index, const size_t nuc, float mass) {
    _seq->update_exp(length()-index-1, complement(nuc), mass); }
  float get_prob(const size_t index, const size_t nuc) const {
    return _seq->get_prob(length()-index-1, complement(nuc)); }
  bool prob() const { return _seq->prob(); }
  void calc_p_vals(std::vector<double>& p_vals) const;
};

#endif

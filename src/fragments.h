/**
 *  fragments.h
 *  express
 *
 *  Created by Adam Roberts on 3/23/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 **/

#ifndef FRAGMENTS_H
#define FRAGMENTS_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <api/BamAlignment.h>
#include "sequence.h"
#include "boost/scoped_ptr.hpp"

typedef size_t TargID;
struct Library;
class Target;
class TargetTable;

/**
 * PairStatus enum.
 * PAIRED denotes that the fragment mapping has both paired end reads.
 * LEFT_ONLY denotes that the single read is not reverse complemented => its
 *           left end is the left fragment end.
 * RIGHT_ONLY denotes that the single read is reverse complemented => its right
 *            end is the right fragment end.
 **/
enum PairStatus { PAIRED, LEFT_ONLY, RIGHT_ONLY };

/**
 * The Indel struct stores the information for a single insertion or deletion.
 *  @author    Adam Roberts
 *  @date      2012
 *  @copyright Artistic License 2.0
 **/
struct Indel {
  /**
   * A public size_t for the position of the Indel in the read. 0-based.
   */
  size_t pos;
  /**
   * A public size_t for the length of the Indel in the read.
   */
  size_t len;
  /**
   * Indel constructor
   */
  Indel(size_t p, size_t l) : pos(p), len(l) {}
};

// DOC
/**
 * A public bool specifying which read comes first in target coordinates.
 * The second read (according to SAM flag) is reverse complemented when
 * true, and the "left" (first according to SAM flag) is reverse complemented
 * when false. In other words, the "left" read is truly left of the "right"
 * read in target coordinate space when true.
 */
struct ReadHit {
  /**
   * A public string for the SAM "Query Template Name" (fragment name)
   */
  std::string name;
  bool first;
  bool reversed;
  /**
   * A public TargID for the target mapped to.
   */
  size_t targ_id;
  /**
   * A public size_t containing the 0-based leftmost coordinate mapped to in the
   * target.
   */
  size_t left;
  /**
   * A public size_t containing the position following the 0-based rightmost
   * coordinate mapped to in the target.
   */
  size_t right;
  /**
   * The read sequence.
   */
  SequenceFwd seq;
  /**
   * A public vector of Indel objects storing all insertions to the reference in
   * the read. Insertions are stored in read order.
   */
  std::vector<Indel> inserts;
  /**
   * A public vector of Indel objects storing all insertions to the reference in
   * the read. Deletions are stored in read order.
   */
  std::vector<Indel> deletes;
  /**
   * A public BamAlignment object storing the raw alignment information from
   * BamTools for the read. Only valid if BAM file is input.
   */
  BamTools::BamAlignment bam;
  mutable std::vector<char> bias_index_cache;
  mutable std::vector<char> mismatch_index_cache;
  /**
   * A public string storing the raw alignment information from for the read.
   * Only valid if SAM file is input.
   */
  std::string sam;
  /**
   * A public int containing the left position for the mate of the read.
   * -1 if single-end fragment. This is temporarily used to help find the mate
   * but is not used after.
   */
  int mate_l;
};

/**
 * The FragHit struct stores the information for a single fragment alignment.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class FragHit {
  /**
   * A public pointer to the target mapped to.
   */
  Target* _target;
  /**
   * Data for the upstream (left) read.
   */
  boost::scoped_ptr<ReadHit> _read_l;
  /**
   * Data for the downstream (right) read.
   * PairStatus is PAIRED.
   */
  boost::scoped_ptr<ReadHit> _read_r;
  /**
   * A public double storing the posterior probability of the mapping after
   * processing. The posteriors of all mappings for a given Fragment should sum
   * to 1.
   */
  double _probability;
  /**
   * A public vector storing pointers to "neighboring" targets. This is being
   * used for an experimental feature and may be removed.
   */
  std::vector<Target*> _neighbors;
  
public:
  FragHit(ReadHit* h) : _target(NULL), _probability(LOG_0){
    if (h->reversed) {
      _read_r.reset(h);
    } else {
      _read_l.reset(h);
    }
  }
  FragHit(ReadHit* l, ReadHit* r) : _read_l(l), _read_r(r) {
    assert(!l->reversed);
    assert(r->reversed);
    assert(l->name == r->name);
    assert(l->targ_id == r->targ_id);
    assert(l->left <= r->left);
    assert(l->first != r->first);
  }
  double probability() const {
    return _probability;
  }
  void probability(double p) {
    _probability = p;
  }
  Target* target() const {
    return _target;
  }
  void target(Target* target) {
    _target = target;
  }
  const std::vector<Target*>* neighbors() const { return &_neighbors; }
  void neighbors(const std::vector<Target*>& neighbors) {
    _neighbors = neighbors;
  }
  TargID target_id() const {
    if (_read_l) {
      return _read_l->targ_id;
    }
    assert(_read_r);
    return _read_r->targ_id;
  }
  size_t left() const {
    if (_read_l) {
      return _read_l->left;
    }
    assert(_read_r);
    return _read_r->left;
  }
  size_t right() const {
    if (_read_r) {
      return _read_r->right;
    }
    assert(_read_l);
    return _read_l->right;
  }
  TargID targ_id() const {
    if (_read_l) {
      return _read_l->targ_id;
    }
    assert(_read_r);
    return _read_r->targ_id;
  }
  /**
   * A member function returning the length of the fragment according to this
   * mapping. Returns 0 if the fragment is single-end.
   * @return Length of fragment mapping.
   */
  size_t length() const {
    if (_read_l && _read_r) {
      return _read_r->right - _read_l->left;
    }
    return 0;
  }
  const ReadHit* left_read() const {
    if (_read_l) {
      return _read_l.get();
    }
    return NULL;
  }
  const ReadHit* right_read() const {
    if (_read_r) {
      return _read_r.get();
    }
    return NULL;
  }
  const ReadHit* first_read() const {
    if (_read_l && _read_l) {
      return _read_l.get();
    }
    assert(_read_r);
    return _read_r.get();
  }
  const ReadHit* second_read() const {
    if (_read_l && !_read_l->first) {
      return _read_l.get();
    } else if (_read_r && !_read_r->first) {
      return _read_r.get();
    } else {
      return NULL;
    }
  }
  ReadHit* left_read() {
    return const_cast<ReadHit*>(const_cast<const FragHit*>(this)->left_read());
  }
  ReadHit* right_read() {
    return const_cast<ReadHit*>(const_cast<const FragHit*>(this)->right_read());
  }
  ReadHit* first_read() {
    return const_cast<ReadHit*>(const_cast<const FragHit*>(this)->first_read());
  }
  ReadHit* second_read() {
    return const_cast<ReadHit*>(const_cast<const FragHit*>(this)->second_read());
  }
  /**
   * A member function returning whether the mapping is PAIRED, LEFT_ONLY, or
   * RIGHT_ONLY, as defined in the PairStatus enum definition.
   * @return The pair status of the mapping
   */
  PairStatus pair_status() const {
    if (_read_l && _read_r) {
      return PAIRED;
    }
    if (_read_l) {
      return LEFT_ONLY;
    }
    return RIGHT_ONLY;
  }
};

/**
 * The Fragment class stores information for all alignments of a single fragment.
 * By design, only paired-end mappings of paired-end reads will be accepted.

 * All mappings of single-end reads will be accepted.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class Fragment {
  /**
   * A private vector of FragHit pointers containing all mappings of the
   * fragment.
   */
  std::vector<FragHit*> _frag_hits;
  /**
   * A private vector of FragHit pointers containing single-end mappings whose
   * pairs have not been found. Temporarily useful when parsing the input to
   * find mates but is not used after.
   */
  std::vector<ReadHit*> _open_mates;
  /**
   * A private string for the SAM "Query Template Name" (fragment name).
   */
  std::string _name;
  /**
   * A private double for the mass of the Fragment as determined by the
   * forgetting factor during processing.
   */
  double _mass;
  /**
   * A private pointer to the global variables associated with the library
   * this fragment is from.
   */
  Library* _lib;
  /**
   * A private method that searches for the mate of the given read mapping.
   * If found, the mates are combined into a single FragHit and added to
   * frag_hits. If not found, the read mapping is added to open_mates.
   */
  void add_open_mate(ReadHit* om);

public:
  /**
   * Fragment Constructor.
   */
  Fragment(Library* lib);
  /**
   * Fragment destructor deletes all FragHit objects pointed to by the Fragment.
   */
  ~Fragment();
  /**
   * Accessor for the global variables associated with the library this fragment
   * is from. Pointer outlives this.
   */
  const Library* lib() { return _lib; }
  /**
   * A member function that adds a new FragHit (single read at this point) to
   * the Fragment. If it is the first FragHit, it sets the Fragment name and is
   * added to _open_mates. If the fragment is not paired, it is added to
   *_frag_hits. Otherwise, add_open_mate is called.
   * @param f the FragHit to be added.
   * @return False iff the read name does not match the Fragment name.
   */
  bool add_map_end(ReadHit* f);
  /**
   * A member function that returns a reference to the "Query Template Name".
   * @return Reference to the SAM "Query Template Name" (fragment name).
   */
  const std::string& name() const { return _name; }
  /**
   * An accessor for the number of valid alignments of the fragment.
   * @return Number of valid alignments for fragment.
   */
  size_t num_hits() const { return _frag_hits.size(); }
  // DOC
  FragHit* operator[](size_t i) const { return _frag_hits[i]; }
  /**
   * Accessor for the FragHit objects associated with the fragment. Returned
   * value does not outlive this.
   * @return Reference to a vector containing pointers to the FragHits.
   */
  const std::vector<FragHit*>& hits() const { return _frag_hits; }
  /**
   * A member function that returns a single FragHit of the fragment sampled at
   * random based on the probabalistic assignments. Returned value does not
   * outlive this.
   * @return A randomly sampled FragHit.
   */
  const FragHit* sample_hit() const;
  /**
   * Mutator for the mass of the fragment according to the forgetting factor.
   * @param m a double representing the value to set to the mass to.
   */
  void mass(double m) { _mass = m; }
  /**
   * An accessor for the mass of the fragment according to the forgetting
   * factor.
   * @return The mass of the fragment.
   */
  double mass() const { return _mass; }
  /**
   * A member function that sorts the FragHits by the TargID of the targets they
   * are aligned to.
   */
  void sort_hits();
};

#endif

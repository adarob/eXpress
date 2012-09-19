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

/**
 * The FragHit struct stores the information for a single fragment alignment.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
struct FragHit {
  /**
   * A public string for the SAM "Query Template Name" (fragment name)
   */
  std::string name;
  /**
   * A public TargID for the target mapped to.
   */
  TargID targ_id;
  /**
   * A public pointer to the target mapped to.
   */
  Target* targ;
  /**
   * The first ("left") read sequence (according to SAM flag).
   */
  SequenceFwd seq_l;
    
  /**
   * The second ("right") read sequence (according to SAM flag).
   */
  SequenceFwd seq_r;
  /**
   * A public size_t containing the 0-based leftmost coordinate mapped to in the
   * target. Valid only iff PairStatus is PAIRED or LEFT_ONLY.
   */
  size_t left;
  /**
   * A public size_t containing the position following the 0-based rightmost
   * coordinate mapped to in the target. Valid iff PairStatus is PAIRED or
   * RIGHT_ONLY.
   */
  size_t right;
  /**
   * A public int containing the left position for the mate of the "left" read.
   * 0 if single-end fragment. This is temporarily used to help find the mate
   * but is not used after.
   */
  int mate_l;
  /**
   * A public bool specifying which read comes first in target coordinates.
   * The "right" (second according to SAM flag) is reverse complemented when
   * true, and the "left" (first according to SAM flag) is reverse complemented
   * when false. In other words, the "left" read is truly left of the "right"
   * read in target coordinate space when true.
   */
  bool left_first;
  /**
   * A public vector of Indel objects storing all insertions to the reference in
   * the "left" read. Insertions are stored in read order.
   */
  std::vector<Indel> inserts_l;
  /**
   * A public vector of Indel objects storing all insertions to the reference in
   * the "left" read. Deletions are stored in read order.
   */
  std::vector<Indel> deletes_l;
  /**
   * A public vector of Indel objects storing all insertions to the reference in
   * the "right" read. Insertions are stored in read order.
   */
  std::vector<Indel> inserts_r;
  /**
   * A public vector of Indel objects storing all deletions to the reference in
   * the "right" read. Deletions are stored in read order.
   */
  std::vector<Indel> deletes_r;
  /**
   * A public double storing the posterior probability of the mapping after
   * processing. The posteriors of all mappings for a given Fragment should sum
   * to 1.
   */
  double probability;
  /**
   * A public BamAlignment object storing the raw alignment information from
   * BamTools for the "left" read. Only valid if BAM file is input.
   */
  BamTools::BamAlignment bam_l;
  /**
   * A public BamAlignment object storing the raw alignment information from
   * BamTools for the "right" read. Only valid if BAM file is input.
   */
  BamTools::BamAlignment bam_r;
  /**
   * A public string storing the raw alignment information from for the "left"
   * read. Only valid if SAM file is input.
   */
  std::string sam_l;
  /**
   * A public string storing the raw alignment information from for the "right"
   * read. Only valid if SAM file is input.
   */
  std::string sam_r;
  /**
   * A public vector storing pointers to "neighboring" targets. This is being
   * used for an experimental feature and may be removed.
   */
  std::vector<Target*> neighbors;
  /**
   * A member function returning the length of the fragment according to this
   * mapping. Note that this result will be invalid if the fragment is
   * single-end.
   * @return Length of fragment mapping.
   */
  size_t length() const {
    assert (right >= left);
    return right - left;
  }
  /**
   * A member function returning whether the mapping is PAIRED, LEFT_ONLY, or
   * RIGHT_ONLY, as defined in the PairStatus enum definition.
   * @return The pair status of the mapping
   */
  PairStatus pair_status() const {
    if (!seq_l.empty() && !seq_r.empty()) {
      return PAIRED;
    }
    if (seq_l.empty()) {
      return RIGHT_ONLY;
    }
    return LEFT_ONLY;
  }
  
  //DOC
  double targ_rho;
  double const_likelihood;
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
  std::vector<FragHit*> _open_mates;
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
  void add_open_mate(FragHit* om);

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
  bool add_map_end(FragHit* f);
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

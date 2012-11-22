/**
 *  directiondetector.h
 *  express
 *
 *  Created by Adam Roberts on 11/21/12.
 *  Copyright 2012 Adam Roberts. All rights reserved.
 *
 */

#ifndef __express__directiondetector__
#define __express__directiondetector__

#include <cstring>

class Fragment;

/**
 * The DirectionDetector class keeps track of the observed fragment directions
 * (forward-reverse or reverse-forward) and whether they are paired or
 * single-end. It can then determine if the numbers in each direction are
 * distributed disproportionately and throw a warning if the proper direction
 * flag has not been specified. Unable to detect when a directional flag has
 * been chosen incorrectly since the incompatible alignments will have been
 * discarded.
 *  @author    Adam Roberts
 *  @date      2012
 *  @copyright Artistic License 2.0
 */
class DirectionDetector {
  /**
   * A private size_t to store the number of paired-end alignments in the
   * forward-reverse direction.
   */
  size_t _num_fr;
  /**
   * A private size_t to store the number of paired-end alignments in the
   * reverse-foward direction.
   */
  size_t _num_rf;
  /**
   * A private size_t to store the number of single-end alignments in the
   * forward direction.
   */
  size_t _num_f;
  /**
   * A private size_t to store the number of single-end alignments in the
   * reverse direction.
   */
  size_t _num_r;
  
public:
  /**
   * The DirectionDetector constructor sets all counts to 0.
   */
  DirectionDetector();
  /**
   * Adds counts for the alignments of the given Fragment.
   * @param f a pointer to the Fragment to count the direction of its
   *        alignments.
   */
  void add_fragment(Fragment* f);
  /**
   * Throws a warning if a disproportionate number of alignments are in one
   * direction and the proper flag has not been specified.
   * @return True iff a warning is thrown.
   */
  bool report_if_improper_direction();
};

#endif /* defined(__express__directiondetector__) */

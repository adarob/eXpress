/**
 *  robertsfilter.h
 *  express
 *
 *  Created by Adam Roberts on 5/1/12.
 *  Copyright (c) 2012 Adam Roberts. All rights reserved.
 */

#ifndef express_robertsfilter_h
#define express_robertsfilter_h

#include <boost/unordered_set.hpp>
#include <queue>
#include <vector>

static size_t DEFAULT_LOC_SIZE = 10000;
static size_t DEFAULT_GLOB_SIZE = 100000;

/**
 * The RobertsFilter class implements a datastructure to test for repeats of
 * a key with high probability, when repeats are most likely to be nearby.
 * Recently observed keys are stored in a local set for a certain number of
 * observations (set by local_size). After this number of observations, it is
 * removed from the local set, and placed in the global set, displacing a random
 * element of this set. To be used when the full set cannot be stored in memory.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class RobertsFilter {
  /**
   * A private queue to store the local keys in FIFO order. Used to know which
   * element to remove from the set next.
   */
  std::queue<std::string> _local_queue;
  /**
   * A private set used to store the local keys. Allows for easy membership
   * testing.
   */
  boost::unordered_set<std::string> _local_set;
  /**
   * A private vector to store the global keys. Used to know which element to
   * randomly replace when the global set is full.
   */
  std::vector<std::string> _global_vector;
  /**
   * A private set used to store the global keys. Allows for easy membership
   * testing.
   */
  boost::unordered_set<std::string> _global_set;
  /**
   * A private size_t specifying the maximum number of keys to store in the
   * local set.
   */
  size_t _local_size;
  /**
   * A private size_t specifying the maximum number of keys to store in the
   * global set.
   */
  size_t _global_size;

 public:
  /**
   * RobertsFilter constructor sets the size of the sets.
   * @param local_size the maximum number of keys to store in the local set.
   * @param global_size the maximum number of keys to store in the global set.
   */
  RobertsFilter(size_t local_size=DEFAULT_LOC_SIZE,
              size_t global_size=DEFAULT_GLOB_SIZE);
  /**
   * A member function that tests for membership of the key in either set. If
   * not found, the key is added to the local set, possibly pushing the oldest
   * key in the local set into the global set, which may displace a global key.
   * @param key the key to be tested for and pushed into the local set.
   * @return True iff the key is in the local or global set.
   */
  bool test_and_push(const std::string& key);
};


#endif

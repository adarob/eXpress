//
//  directiondetector.h
//  express
//
//  Created by Adam Roberts on 11/21/12.
//
//

#ifndef __express__directiondetector__
#define __express__directiondetector__

#include <cstring>

class Fragment;

class DirectionDetector {
  size_t _num_fr;
  size_t _num_rf;
  size_t _num_f;
  size_t _num_r;
  
public:
  DirectionDetector();
  void add_fragment(Fragment* f);
  bool report_if_improper_direction();
};

#endif /* defined(__express__directiondetector__) */

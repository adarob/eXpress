//
//  directiondetector.cpp
//  express
//
//  Created by Adam Roberts on 11/21/12.
//
//

#include "directiondetector.h"

#include <iostream>
#include "fragments.h"
#include "main.h"

using namespace std;

DirectionDetector::DirectionDetector()
: _num_fr(0), _num_rf(0), _num_f(0), _num_r(0) {}

void DirectionDetector::add_fragment(Fragment* f) {
  foreach (FragHit* h, f->hits()) {
    switch(h->pair_status()) {
      case PAIRED: {
        if (h->left_read() == h->first_read()) {
          _num_fr++;
        } else {
          assert(h->right_read() == h->first_read());
          _num_rf++;
        }
        break;
      }
      case LEFT_ONLY: {
        _num_f++;
        break;
      }
      case RIGHT_ONLY: {
        _num_r++;
        break;
      }
    }
  }
}

bool DirectionDetector::report_if_improper_direction() {
  size_t num_paired = _num_fr + _num_rf;
  size_t num_single = _num_f + _num_r;
  if (num_paired + num_single == 0) {
    return false;
  }
  if (num_paired == 0) {
    // Single-end case
    double max_dir = max(_num_f, _num_r);
    double min_dir = min(_num_f, _num_r);
    if (min_dir < max_dir / 2) {
      if (_num_f > _num_r && direction != F) {
        cerr << "WARNING: The observed alignments appear disporportionately on "
             << "the forward strand (" << _num_f << " vs. " << _num_r << "). "
             << "If your library is strand-specific and single-end, you should "
             << "use the --f-stranded option to avoid incorrect results.\n";
        return true;
      } else if (_num_f < _num_r && direction != R) {
        cerr << "WARNING: The observed alignments appear disporportionately on "
             << "the reverse strand (" << _num_r << " vs. " << _num_f << "). "
             << "If your library is strand-specific and single-end, you should "
             << "use the --r-stranded option to avoid incorrect results.\n";
        return true;
      }
    }
  } else {
    // Paired-end case
    size_t fr = _num_f + _num_fr;
    size_t rf = _num_r + _num_rf;
    double max_dir = max(fr, rf);
    double min_dir = min(fr, rf);
    if (min_dir < max_dir / 2) {
      if (fr > rf && direction != FR) {
        cerr << "WARNING: The observed alignments appear disporportionately in "
        << "the forward-reverse order (" << fr << " vs. " << rf << "). "
        << "If your library is strand-specific, you should use the "
        << "--fr-stranded option to avoid incorrect results.\n";
        return true;
      } else if (rf < fr && direction != RF) {
        cerr << "WARNING: The observed alignments appear disporportionately in "
        << "the reverse-forward order (" << rf << " vs. " << fr << "). "
        << "If your library is strand-specific, you should use the "
        << "--rf-stranded option to avoid incorrect results.\n";
        return true;
      }
    }
  }
  return false;
}

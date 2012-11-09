#ifndef FREQUENCYMATRIX_H
#define FREQUENCYMATRIX_H

/**
 *  frequencymatrix.h
 *  express
 *
 *  Created by Adam Roberts on 4/23/11.
 *  Copyright 2011 Adam Roberts. All rights reserved.
 **/

#include <cassert>
#include <vector>
#include "main.h"

/**
 * The FrequencyMatrix class keeps track of the frequency parameters in order to
 * allow for constant-time probability look-ups and updates. The table is 2D
 * to allow multiple distributions to be stored in one FrequencyMatrix. The
 * first dimension (rows) are the different distributions. Values are stored in
 * log space by default.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
template <class T>
class FrequencyMatrix {
  /**
   * A private vector to store the matrix frequencies (logged) in row-major
   * format.
   */
  std::vector<T> _array;
  /**
   * A private vector to store the (logged) row sums for the matrix. This is the
   * normalizer for the distribution.
   */
  std::vector<T> _rowsums;
  /**
   * A private size_t for the number of rows (distributions).
   */
  size_t _M;
  /**
   * A private size_t for the number of columns.
   */
  size_t _N;
  /**
   * A private bool that specifies if the table values are logged.
   */
  bool _logged;
  /**
   * A private bool that specifies whether or not the values are fixed (locked)
   * and normalized.
   */
  bool _fixed;

public:
  /**
   * Dummy FrequencyMatrix Constructor.
   */
  FrequencyMatrix(){};
  /**
   * FrequencyMatrix constructor initializes the matrix values to the given
   * pseudo-counts.
   * @param m a size_t specifying the number of distributions (rows).
   * @param n a size_t specifying the number of values in each distribution
            (columns).
   * @param alpha the intial psuedo-counts (un-logged).
   * @param logged bool that specifies if the table is to be stored logged.
   */
  FrequencyMatrix(size_t m, size_t n, T alpha, bool logged = true);
  /**
   * An accessor for the frequency at a given position in the matrix (logged
   * if table is logged).
   * @param i the distribution (row).
   * @param j the value (column).
   * @param normalized a bool specifying whether or not the frequency should be
   *        normalized.
   * @return The frequency of the given value in the given distribution
   *         (logged if table is logged).
   */
  T operator()(size_t i, size_t j, bool normalized=true) const;
  /**
   * An accessor for the frequency at a given position in the flattened matrix
   * (logged if table is logged).
   * @param k the array position.
   * @param normalized a bool specifying whether or not the frequency should be
   *        normalized.
   * @return The frequency at the given position in the flattened matrix (logged
   *         if table is logged).
   */
  T operator()(size_t k, bool normalized=true) const;
  /**
   * A member function to increase the mass of a given position in the matrix.
   * @param i the distribution (row).
   * @param j the value (column).
   * @param incr_amt the amount to increase the mass by (logged if table is
   *        logged).
   */
  void increment(size_t i, size_t j, T incr_amt);
  /**
   * A member function to increase the mass of a given position in the flattened
   * matrix (logged if table is logged). Does nothing if _fixed is true.
   * @param k the array position.
   * @param incr_amt the amount to increase the mass by (logged if table is
   *        logged).
   */
  void increment(size_t k, T incr_amt);
  /**
   * A member function to decrease the mass of a given position in the matrix.
   * @param i the distribution (row).
   * @param j the value (column).
   * @param incr_amt the amount to decrease the mass by (logged if table is
   *        logged).
   */
  void decrement(size_t i, size_t j, T decr_amt);
  /**
   * A member function to decrease the mass of a given position in the flattened
   * matrix (logged if table is logged). Does nothing if _fixed is true.
   * @param k the array position.
   * @param incr_amt the amount to decrease the mass by (logged if table is
   *        logged).
   */
  void decrement(size_t k, T decr_amt);
  /**
   * An accessor for the row sum (normalizer), (logged if table is logged).
   * @param i the distribution (row).
   * @return The sum (normalizer) for the given distribution (logged if table is
   *         (logged).
   */
  T sum(size_t i) const { return _rowsums[i]; }
  /**
   * A member function that finds and returns the argmax (index of mode) of the
   * given distribution.
   * @param i the distribution (row).
   * @return The argmax of the distribution.
   */
  size_t argmax(size_t i) const;
  /**
   * A member function that converts the table between log-space and non-log
   * space. Does nothing if _fixed is true.
   * @param logged bool specifying if the table should be converted to logged or
   *        non-logged space.
   */
  void set_logged(bool logged);
  /**
   * A member function that normalizes and "locks" the matrix values so that
   * no changes can be made. This allows for faster future lookups, but is
   * irrevocable.
   */
  void fix();
  /**
   * An accessor for the value of _fixed, which specifies whether or not the
   * matrix has been fixed (irrevocable).
   */
  bool is_fixed() const { return _fixed; }
};

template <class T>
FrequencyMatrix<T>::FrequencyMatrix(size_t m, size_t n, T alpha, bool logged)
    : _array(m*n, logged ? log(alpha):alpha),
      _rowsums(m, logged ? log(n*alpha):n*alpha),
      _M(m),
      _N(n),
      _logged(logged),
      _fixed(false){
}

template <class T>
T FrequencyMatrix<T>::operator()(size_t i, size_t j, bool normalized) const {
  assert(i*_N+j < _M*_N);
  if (_fixed || !normalized) {
      return _array[i*_N+j];
  }
  if (_logged) {
    return _array[i*_N+j]-_rowsums[i];
  } else {
    return _array[i*_N+j]/_rowsums[i];
  }
}

template <class T>
T FrequencyMatrix<T>::operator()(size_t k, bool normalized) const {
  return operator()(0, k, normalized);
}

template <class T>
void FrequencyMatrix<T>::increment(size_t i, size_t j, T incr_amt) {
  if (_fixed) {
    return;
  }

  size_t k = i*_N+j;
  assert(k < _M*_N);
  if (_logged) {
    _array[k] = log_add(_array[k], incr_amt);
    _rowsums[i] = log_add(_rowsums[i], incr_amt);
  } else {
    _array[k] += incr_amt;
    _rowsums[i] += incr_amt;
  }
  assert(!std::isnan(_rowsums[i]));
}

template <class T>
void FrequencyMatrix<T>::decrement(size_t i, size_t j, T decr_amt) {
  if (_fixed) {
    return;
  }
  
  size_t k = i*_N+j;
  assert(k < _M*_N);
  if (_logged) {
    _array[k] = log_sub(_array[k], decr_amt);
    _rowsums[i] = log_sub(_rowsums[i], decr_amt);
  } else {
    _array[k] -= decr_amt;
    _rowsums[i] -= decr_amt;
  }
  assert(!std::isnan(_rowsums[i]));
}

template <class T>
void FrequencyMatrix<T>::increment(size_t k, T incr_amt) {
  increment(0, k, incr_amt);
}

template <class T>
void FrequencyMatrix<T>::decrement(size_t k, T decr_amt) {
  decrement(0, k, decr_amt);
}

template <class T>
void FrequencyMatrix<T>::set_logged(bool logged) {
  if (logged == _logged || _fixed) {
    return;
  }
  if (logged) {
    for (size_t i = 0; i < _M*_N; ++i) {
      _array[i] = log(_array[i]);
    }
    for(size_t i = 0; i < _M; ++i) {
      _rowsums[i] = log(_rowsums[i]);
    }
  } else {
    for (size_t i = 0; i < _M*_N; ++i) {
       _array[i] = sexp(_array[i]);
    }
    for (size_t i = 0; i < _M; ++i) {
      _rowsums[i] = sexp(_rowsums[i]);
    }
  }
  _logged = logged;
}

template <class T>
size_t FrequencyMatrix<T>::argmax(size_t i) const {
  size_t k = i*_N;
  size_t arg = 0;
  T val = _array[k];
  for (size_t j = 1; j < _N; j++) {
    if (_array[k+j] > val) {
      val = _array[k+j];
      arg = j;
    }
  }
  return arg;
}

template <class T>
void FrequencyMatrix<T>::fix() {
  if (_fixed) {
    return;
  }
  for (size_t i = 0; i < _M; ++i) {
    for (size_t j = 0; j < _N; ++j) {
      _array[i*_N+j] = operator()(i,j);
    }
    _rowsums[i] = (_logged) ? 0 : 1;
  }
  _fixed = true;
}

#endif

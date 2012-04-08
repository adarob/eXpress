#ifndef FREQUENCYMATRIX_H
#define FREQUENCYMATRIX_H

//
//  frequencymatrix.h
//  express
//
//  Created by Adam Roberts on 4/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//
#include <vector>
#include <cassert>
#include "main.h"

/**
 * The FrequencyMatrix class keeps track of the frequency parameters
 * in order to allow for constant-time probability look-ups and updates.
 * The table is rectangular to allow for multiple distributions to be
 * stored in one FrequencyMatrix.  Rows are distributions. Values
 * are in log space by default.
 * @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
template <class T>
class FrequencyMatrix
{
   /**
    * a private vector of doubles to store the matrix frequencies (logged) in row-major format
    */
    std::vector<T> _array;
    
   /**
    * a private vector of doubles to store the (logged) row sums for the matrix
    */
    std::vector<T> _rowsums;
    
    /**
     * a private size_t for the number of rows
     */
    size_t _M;
    
    /**
     * a private size_t for the number of columns
     */
    size_t _N;
    
    /**
     * a private bool that specifies if the table is in log space
     */
    bool _logged;
    
public:
    
    /**
     * dummy constructor
     */    
    FrequencyMatrix(){};
    
    /**
     * FrequencyMatrix constructor initializes the matrix based on the log of the given pseudo-counts
     * @param m a size_t specifying the number of distributions (rows)
     * @param n a size_t specifying the number of values in each distribution (columns)
     * @param alpha a double specifying the intial psuedo-counts (un-logged)
     * @param logged bool that specifies if the table is in log space
     */   
    FrequencyMatrix(size_t m, size_t n, T alpha, bool logged = true);
   
    /**
     * a member function to extract the probability of a given position in the matrix (logged if table is logged)
     * @param i the distribution (row)
     * @param j the value (column)
     * @return a double specifying the probability of the given value in the given distribution (logged if table is logged)
     */  
    T operator()(size_t i, size_t j) const;
    
    /**
     * a member function to extract the probability of a given position in the flattened matrix (logged if table is logged)
     * @param k the array position
     * @return a double specifying the probability of the given position in the flattened matrix (logged if table is logged)
     */  
    T operator()(size_t k) const;
    
    /**
     * a member function to increase the mass of a given position in the matrix
     * @param i the distribution (row)
     * @param j the value (column)
     * @param incr_amt the amount to increase the mass by (logged if table is logged)
     */ 
    void increment(size_t i, size_t j, T incr_amt);
    
    /**
     * a member function to increase the mass of a given position in the flattened matrix (logged if table is logged)
     * @param k the array position
     * @param incr_amt the amount to increase the mass by (logged if table is logged)
     */ 
    void increment(size_t k, T incr_amt);
    
    /**
     * a member function that returns the raw value stored at a given position of the flattened matrix
     * @param k the array position
     * @return a double specifying the raw value stored at the given position of the flattened matrix
     */ 
    double arr(size_t k) const { return _array[k]; }
    
    /**
     * a member function that returns the raw row sum
     * @param i the distribution (row)
     * @return a double specifying the raw row sum for the given distribution
     */ 
    double row(size_t i) const { return _rowsums[i]; } 
    
    /**
     * a member function that converts the table between log-space and non-log space
     * @param logged bool specifying if the table should be converted to logged or non-logged space
     */ 
    void set_logged(bool logged);
};


template <class T>
FrequencyMatrix<T>::FrequencyMatrix(size_t m, size_t n, T alpha, bool logged)
: _array(m*n, logged ? log(alpha):alpha),
_rowsums(m, logged ? log(n*alpha):n*alpha),
_M(m),
_N(n),
_logged(logged)
{}

template <class T>
T FrequencyMatrix<T>::operator()(size_t i, size_t j) const
{
    assert(i*_N+j < _M*_N);
    if (_logged)
        return _array[i*_N+j]-_rowsums[i];
    else
        return _array[i*_N+j]/_rowsums[i];
}

template <class T>
T FrequencyMatrix<T>::operator()(size_t k) const
{
    return operator()(0, k);
}

template <class T>
void FrequencyMatrix<T>::increment(size_t i, size_t j, T incr_amt)
{
    size_t k = i*_N+j;
    assert(k < _M*_N);
    if (_logged)
    {
        _array[k] = log_sum(_array[k], incr_amt);
        _rowsums[i] = log_sum(_rowsums[i], incr_amt);
    }
    else
    {
        _array[k] += incr_amt;
        _rowsums[i] += incr_amt;
    }
//    assert(!isnan(_rowsums[i]) && !isinf(_rowsums[i]));
}

template <class T>
void FrequencyMatrix<T>::increment(size_t k, T incr_amt)
{
    increment(0, k, incr_amt);
}

template <class T>
void FrequencyMatrix<T>::set_logged(bool logged)
{
    if (logged == _logged)
        return;
    
    if (logged)
    {
        for(size_t i = 0; i < _M*_N; ++i)
        {
            _array[i] = log(_array[i]);
        }
        for(size_t i = 0; i < _M; ++i)
        {
            _rowsums[i] = log(_rowsums[i]);
        }
    }
    else
    {
        for(size_t i = 0; i < _M*_N; ++i)
        {
            _array[i] = sexp(_array[i]);
        }
        for(size_t i = 0; i < _M; ++i)
        {
            _rowsums[i] = sexp(_rowsums[i]);
        }
    }
    
    _logged = logged;    
}

#endif
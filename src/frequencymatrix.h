#ifndef FREQUENCYMATRIX_H
#define FREQUENCYMATRIX_H

//
//  frequencymatrix.h
//  express
//
//  Created by Adam Roberts on 4/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//
#include "boost/multi_array.hpp"
#include <vector>

/**
 * The FrequencyMatrix class keeps track of the frequency parameters
 * in order to allow for constant-time probability look-ups and updates.
 * The table is rectangular to allow for multiple distributions to be
 * stored in one FrequencyMatrix.  Rows are distributions. All values 
 * are stored and returned in log space.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class FrequencyMatrix
{
   /**
    * a private vector of doubles to store the matrix frequencies (logged) in row-major format
    */
    std::vector<double> _array;
    
   /**
    * a private vector of doubles to store the (logged) row sums for the matrix
    */
    std::vector<double> _rowsums;
    
    /**
     * a private size_t for the number of rows
     */
    size_t _M;
    
    /**
     * a private size_t for the number of columns
     */
    size_t _N;
    
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
     */   
    FrequencyMatrix(size_t m, size_t n, double alpha);
   
    /**
     * a member function to extract the logged probability of a given position in the matrix
     * @param i the distribution (row)
     * @param j the value (column)
     * @return a double specifying the logged probability of the given value in the given distribution
     */  
    double operator()(size_t i, size_t j) const;
    
    /**
     * a member function to extract the logged probability of a given position in the flattened matrix
     * @param k the array position
     * @return a double specifying the logged probability of the given position in the flattened matrix
     */  
    double operator()(size_t k) const;
    
    /**
     * a member function to increase the mass of a given position in the matrix
     * @param i the distribution (row)
     * @param j the value (column)
     * @param incr_amt the logged amount to increase the mass by
     */ 
    void increment(size_t i, size_t j, double incr_amt);
    
    /**
     * a member function to increase the mass of a given position in the flattened matrix
     * @param k the array position
     * @param incr_amt the logged amount to increase the mass by
     */ 
    void increment(size_t k, double incr_amt);
    
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
};

#endif
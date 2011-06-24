//
//  fld.h
//  expressionline
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#ifndef FLD_H
#define FLD_H

#include <vector>

class FLD
{
    /**
     * a private vector that stores the observed counts by length
     */
    std::vector<long double> _hist;
    
    /**
     * a private long double that stores the total observed counts
     */
    long double _num_obs;
    
    /**
     * a private long double that stores the sum of the observations for quick mean calculations
     */
    long double _sum;
    
    /**
     * a private int that stores the minimum observed length
     */
    size_t _min;
    
public:
    /**
     * FLD Constructor
     * @param alpha an integer that sets the initial uniform pseudo-counts
     * @param max_val an integer that sets the maximum allowable FragMap length
     */
    FLD(long double alpha, size_t max_val);
    
    /**
     * a member function that returns the maximum allowed FragMap length
     * @return max allowed FragMap length
     */ 
    size_t max_val() const;
    
    /**
     * a member function that returns the mean FragMap length
     * @return mean observed FragMap length
     */ 
    long double mean() const;
    
    /**
     * a member function that updates the distribution based on a new FragMap observation
     * @param length an integer for the observed FragMap length
     * @param mass a long double for the mass of the observed FragMap
     */
    void add_val(size_t len, long double mass);
    
    /**
     * a member function that returns the probability of a given FragMap length
     * @param length an integer for the FragMap length to return the probability of
     * @return probability of observing the given FragMap length
     */ 
    long double pdf(size_t len) const;
    
    /**
     * a member function that returns the probability of any length less than 
     * or equal to the given FragMap length
     * @param length an integer for the length to return probability of a FragMap being less than or equal to
     * @return probabilty of observing a FragMap length less than or equal to the given length
     */
    long double cdf(size_t len) const;
    
    /**
     * a member function that returns the probability of the given length conditioned on the max value
     * @param length an integer for the FragMap length to return renormalized probability of
     * @param max an integer specifying the max nonzero probability length to renormalize the distribution
     * @return probabilty of observing the given FragMap length conditioning on the maximum value
     */
    long double npdf(size_t len, size_t max) const;
    
    /**
     * a member function that returns the number of observed fragmaps (including pseudo-counts)
     * @return number of observed fragments
     */ 
    long double num_obs() const;
    
    void output(std::string path) const;

};

#endif
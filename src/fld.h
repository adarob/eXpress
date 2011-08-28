//
//  fld.h
//  express
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#ifndef FLD_H
#define FLD_H

#include <vector>
#include <string>

/**
 * The FLD class keeps track of the observed fragment length distribution.  It starts with
 * a Gaussian prior with parameters specified by the arguments.  A small "Gaussian" kernel is added
 * for each observation.  All mass values and probabilities are stored and returned in log space 
 * (except in to_string).
 */
class FLD
{
    /**
     * a private vector that stores the observed counts (logged) by length
     */
    std::vector<double> _hist;
    
    /**
     * a private double that stores the total observed counts (logged)
     */
    double _tot_mass;
    
    /**
     * a private double that stores the (logged) sum of the observed lengths for quick mean calculations
     */
    double _sum;
    
    /**
     * a private int that stores the minimum observed length
     */
    size_t _min;
    
public:
    /**
     * FLD Constructor
     * @param alpha double that sets the average pseudo-counts (logged)
     * @param max_val an integer that sets the maximum allowable FragMap length
     * @param mean a size_t for the mean of the prior gaussian dist
     * @param std_dev a size_t for the std dev of the prior gaussian dist
     */
    FLD(double alpha, size_t max_val, size_t mean, size_t std_dev);
    
    /**
     * a member function that returns the maximum allowed FragMap length
     * @return max allowed FragMap length
     */ 
    size_t max_val() const;
    
    /**
     * a member function that returns the mean FragMap length
     * @return mean observed FragMap length
     */ 
    double mean() const;
    
    /**
     * a member function that updates the distribution based on a new FragMap observation
     * @param len an integer for the observed FragMap length
     * @param mass a double for the mass (logged) of the observed FragMap
     */
    void add_val(size_t len, double mass);
    
    /**
     * a member function that returns the (logged) probability of a given FragMap length
     * @param len an integer for the FragMap length to return the probability of
     * @return (logged) probability of observing the given FragMap length
     */ 
    double pdf(size_t len) const;
    
    /**
     * a member function that returns the (logged) number of observed fragmaps (including pseudo-counts)
     * @return number of observed fragments
     */ 
    double tot_mass() const;
    
    /**
     * a member function that returns a string containing the current distribution
     * @return space-separated string of probabilities ordered from length 0 to max_val (non-logged)
     */ 
    std::string to_string() const;

};

#endif

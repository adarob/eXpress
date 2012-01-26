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
     * @param max_val an integer that sets the maximum allowable FragHit length
     * @param mean a size_t for the mean of the prior gaussian dist
     * @param std_dev a size_t for the std dev of the prior gaussian dist
     */
    FLD(double alpha, size_t max_val, size_t mean, size_t std_dev);
    
    /**
     * a member function that returns the maximum allowed FragHit length
     * @return max allowed FragHit length
     */ 
    size_t max_val() const;
    
    /**
     * a member function that returns the mean FragHit length
     * @return mean observed FragHit length
     */ 
    double mean() const;
    
    /**
     * a member function that updates the distribution based on a new FragHit observation
     * @param len an integer for the observed FragHit length
     * @param mass a double for the mass (logged) of the observed FragHit
     */
    void add_val(size_t len, double mass);
    
    /**
     * a member function that returns the (logged) probability of a given FragHit length
     * @param len an integer for the FragHit length to return the probability of
     * @return (logged) probability of observing the given FragHit length
     */ 
    double pdf(size_t len) const;
    
    /**
     * a member function that returns the (logged) number of observed FragHits (including pseudo-counts)
     * @return number of observed fragments
     */ 
    double tot_mass() const;
    
    /**
     * a member function that returns a string containing the current distribution
     * @return space-separated string of probabilities ordered from length 0 to max_val (non-logged)
     */ 
    std::string to_string() const;

    /**
     * a member function that appends the FLD parameters to the end of the given file
     * @param outfile the file to append to
     */ 
    void append_output(std::ofstream& outfile) const;
    
};

#endif

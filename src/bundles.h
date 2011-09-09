//
//  bundles.h
//  express
//
//  Created by Adam Roberts on 9/6/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#ifndef express_bundles_h
#define express_bundles_h

#include <vector>
#include <boost/unordered_set.hpp>


class Transcript;
class FragMassTable;

/** 
 * The Bundle class keeps track of a group of transcripts that have shared ambiguous (multi-mapped)
 * reads. Besides storing the transcript, it keeps track of the number of observed fragments,
 * the total fragment mass, and the next fragment mass (which it also updates).
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class Bundle
{
    /**
     * a private vector that stores pointers to all transcripts in the bundle
     */
    std::vector<Transcript*> _transcripts;     
    
    /**
     * a private double that stores the total mass of all fragments mapped to the bundle (logged)
     */
    double _mass;
    
    /**
     * a private size_t that stores the total number of observed fragments mapped to the bundle
     */
    size_t _counts;
    
    /**
     * a private double that stores the total number of pseudo counts for transcripts in the bundle (logged)
     */
    double _pseudo_mass;
    
    /**
     * a private size_t that acts as proxy for observed counts, allowing stored mass and cumulative mass values to be used in the FragmentMassTable
     */
    size_t _n;
    
    /**
     * a private double that stores the mass of the next observed fragment in the bundle
     */
    double _next_frag_mass;
    
    /**
     * a private pointer to the FragMassTable which calculates the fragment mass values
     */
    FragMassTable* _fmt;
    
public:
    
    /**
     * Bundle Constructor.
     * @param trans a pointer to the initial Transcript object in the bundle
     * @param fmt a pointer to the (global) FragMassTable
     */
    Bundle(Transcript* trans, FragMassTable* fmt);
    
    /**
     * a member function that returns the mass of the next observed fragment and has the FragMassTable calculate the subsequent value
     * @return the mass of the next observed fragment for the bundle (logged)
     */
    double next_frag_mass();
    
    /**
     * a member function that sets the value of the mass for the next observed fragment in the bundle, presumably after a merge 
     * @param new value for the mass of the next observed fragment for the bundle (logged)
     */
    void next_frag_mass(double m) { _next_frag_mass = m; }
    
    /**
     * a member function that sets the value of the observed fragment count proxy, presumably after a merge 
     * @param new value for the observed fragment count proxy
     */
    void n(size_t n) { _n = n; }
    
    /**
     * a member function that sets the total fragment mass of the bundle, presumably after a merge
     * @ the new fragment mass of the bundle (logged)
     */
    void mass(double m) { _mass = m; }
    
    
    /**
     * a member function that renormalizes the mass of all transcripts in the bundle by multiplying (adding in log space) the normalizing constant
     * @param norm_const the constant to multiply (add in log space) the transcript masses by (logged)
     */
    void renormalize_transcripts(double norm_const);
    
    /**
     * a member function that increases the total bundle mass by a given amount
     * @param mass the amount to increase the bundle mass by (logged)
     */
    void add_mass(double mass);
    
    /**
     * a member function that increases the total psuedo mass by a given amount
     * @param mass the amount to increase the bundle pseudo mass by (logged)
     */
    void add_pseudo_mass(double mass);
    
    /**
     * a member function that increases the total bundle observed fragment counts by a given amount
     * @param incr_amt the amount to increase the counts by
     */
    void incr_counts(size_t incr_amt=1);
    
    /**
     * a member function that returns the number of transcripts in the bundle
     * @return the number of transcripts in the bundle
     */
    size_t size() const { return _transcripts.size(); }
    
    /**
     * a member function that returns a reference to the vector of pointers to transcripts in the bundle
     * @return reference to the vector pointing to bundle transcripts
     */
    std::vector<Transcript*>& transcripts() { return _transcripts; }
    
    /**
     * a member function that returns the total fragment mass of the bundle
     * @return the total fragment mass of the bundle (logged)
     */
    double mass() const { return _mass; }
    
    /**
     * a member function that returns the total pseudo-mass for transcripts in the bundle (logged)
     * @return the total pseudo-mass of the bundle (logged)
     */
    double pseudo_mass() const { return _mass; }
    
    /**
     * a member function that returns the total number of observed fragments mapped to transcripts in the bundle
     * @return the total number of fragments mapped to transcripts in the bundle
     */
    size_t counts() const { return _counts; }
    
    /**
     * a member function that returns the total number of observed fragments mapped to transcripts in the bundle
     * plus the psuedo-counts weighted by their proportion of the current mass
     * @return the total number of fragments mapped to bundle plus weighted psuedo-counts
     */
    double counts_plus_pseudo() const;
};


typedef boost::unordered_set<Bundle*> BundleSet;

/** 
 * The BundleTable class keeps track of the Bundle objects for a given run.  It has the ability to create, delete,
 * and merge bundles. It also keeps track of the transcript covariances, since these are related to bundles in that
 * all covariances outside of a bundle are nonzero.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class BundleTable
{
    /**
     * a private unordered_set to store all of the bundles
     */
    BundleSet _bundles;
    
    /**
     * a private pointer to the (global) FragMassTable to pass to Bundle objects and help renormalize during a merge
     */
    FragMassTable* _fmt;

public:
    
    /**
     * BundleTable constructor.
     * @param fmt a pointer to the (global) FragMassTable
     */
    BundleTable(FragMassTable* fmt);
    
    /**
     * BundleTable deconstructor.  Deletes all Bundle objects.
     */
    ~BundleTable();
    
    /**
     * a member function that returns the set of current Bundle objects
     * @return a reference to the unordered_set containing all current Bundle objects
     */
    const BundleSet& bundles() const { return _bundles; }
    
    /**
     * a member function that returns the size of the set of current Bundle objects,
     * which is the current number of bundles
     * @return the current number of bundles
     */    
    size_t size() const { return _bundles.size(); }
    
    /**
     * a member function that creates a new bundle, initially with only the single given Transcript
     * @param trans a pointer to the only Transcript initially contained in the Bundle
     * @return a pointer to the new Bundle object
     */ 
    Bundle* create_bundle(Transcript* trans);
    
    /**
     * a member function that merges two Bundle objects into one
     * a new total mass for the merged Bundle is calculated based on the combined number of observed fragments
     * the masses of the member Transcripts are renormalized to proportionally distribute their masses in a way that sums to this new total mass
     * the Transcripts all shifted to one of the bundles (the larger one) and the other is deleted
     * @param b1 a pointer to one of the Bundle objects to merge
     * @param b2 a pointer to the other Bundle object to merge
     * @return a pointer to the merged Bundle object
     */ 
    Bundle* merge(Bundle* b1, Bundle* b2);
};


#endif

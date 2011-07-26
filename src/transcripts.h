//
//  transcripts.h
//  expressionline
//
//  Created by Adam Roberts on 3/20/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#ifndef TRANSCRIPTS_H
#define TRANSCRIPTS_H

#include <string>
#include <map>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include <vector>
#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include "main.h"


class FLD;
class FragMap;
class BiasBoss;
class MismatchTable;

typedef size_t TransID;

static boost::hash<std::string> hash_trans_name;

static size_t hash_trans_pair(TransID trans1, TransID trans2)
{
    size_t seed = std::min(trans1, trans2);
    boost::hash_combine(seed, std::max(trans1, trans2));
    return seed;
}

/**
 * Transcript class.  This class is used to store objects for the transcripts
 * being mapped to.  Besides storing basic information about the object
 * (id, length), it also stores the current expected the number of 
 * transcripts mapping to the object.  To help with updating this number,
 * it will also return the likelihood that a given fragment originated from it.
 */
class Transcript
{
    /**
     * a private string that stores the transcript name
     */
    std::string _name;
    
    TransID _id;
    
    std::string _seq;
    
    /**
     * a private unsigned integer that stores the transcript length
     */
    size_t _len;
    
    /**
     * a private double that stores the expected counts
     */
    double _counts;
    
    double _var;
    
    mutable boost::mutex _bias_lock;
    std::vector<double> _start_bias;
    std::vector<double> _end_bias;
    double _avg_bias;
    
    /**
     * a pointer to the Fragment Length Distribution (FLD) object
     */
    const FLD* _fld;
    
    const BiasBoss* _bias_table;
    
    const MismatchTable* _mismatch_table;
        
    double total_bias_for_length(size_t l) const;
    
public:
    
    /**
     * Transcript Constructor
     * @param id a string that stores the transcript name/id
     * @param seq a string that stores the transcript sequence
     * @param alpha a double that specifies the intial pseudo-counts
     * @param fld a pointer to the Fragment Length Distribution (FLD) object
     */
    Transcript(const std::string& name, const std::string& seq, double alpha, const FLD* fld, const BiasBoss* bias_table, const MismatchTable* mismatch_table);
    
    /**
     * a member function that returns the transcript name
     * @return transcript id
     */
    const std::string& name() const { return _name; }
    
    /**
     * a member function that returns the transcript id
     * @return transcript name
     */
    TransID id() const { return _id; }
    
    const std::string& seq() const { return _seq; }
    /**
     * a member function that returns the transcript length
     * @return transcript length
     */
    size_t length() const { return _len; }
    
    /**
     * a member function that returns the current expected fragment count
     * @return expected fragment count
     */
    double frag_count() const { return _counts; }
    
    /**
     * a member function that increases the expected fragment counts by a given mass
     * @param mass a double to increase the expected fragment counts by
     * @param the probabilistic assignment of the 
     */
    void add_mass(double p, double mass) 
    { 
        assert(sexp(p) != 0);
        _counts = log_sum(_counts, p+mass);
        if (p != 0.0)
            _var = log_sum(_var, mass + p + log(1-sexp(p)));
    }  
    
    void set_counts(double counts)
    {
        _counts = counts;
    }
    
    /**
     * a member function that returns (a value proportional to) the likelihood the given fragment
     * originated from this transcript
     * @param frag a FragMap to return the likelihood of being originated from this transcript
     */
    double log_likelihood(const FragMap& frag) const;

    double effective_length() const;
    
    void update_transcript_bias();
    
    const FLD* fld() { return _fld; }
};

typedef boost::unordered_map<TransID, Transcript*> TransMap;
typedef boost::unordered_map<size_t, double> TransPairMap;

/**
 * TranscriptTable class.  This class is used to keep track of the Transcript objects for a run.
 * The constructor parses a fasta file to generate the Transcript objects and store them in a map
 * that allows them to be looked up based on their string id.
 */
class TranscriptTable
{
    /**
     * a map to look up pointers to Transcript objects by their string id
     */
    TransMap _trans_map;
    TransPairMap _covar_map;
    double _alpha;
    
    /**
     * a private function that adds a transcript pointer to the table
     * @param trans a pointer to the transcript to be added
     */
    void add_trans(Transcript* trans);  
public:
    /**
     * TranscriptTable Constructor
     * @param trans_fasta_file a string storing the path to the fasta file from which to load transcripts
     * @param alpha a double that specifies the intial pseudo-counts for the transcripts
     * @param fld a pointer to the Fragment Length Distribution (FLD) object
     */
    TranscriptTable(const std::string& trans_fasta_file, double alpha, const FLD* fld, BiasBoss* bias_table, const MismatchTable* mismatch_table);
    
    /**
     * TranscriptTable Destructor
     * deletes all of the transcript objects in the table
     */
    ~TranscriptTable();
    
    /**
     * a member function that returns a pointer to the transcript with the given id 
     * @param trans_name id of the transcript wanted
     * @return pointer to the transcript wit hthe given id
     */
    Transcript* get_trans(TransID id);
    
    /**
     * a member function that returns the number of transcripts in the table
     * @return number of transcripts in the table
     */
    size_t size() const { return _trans_map.size(); }
    
    void update_covar(TransID trans1, TransID trans2, double covar);
    size_t covar_size() const { return _covar_map.size(); }
    
    void threaded_bias_update();
    
    void output_header(std::ofstream& runexpr_file);
    void output_current(std::ofstream& runexpr_file);

    void output_expression(std::string output_dir);
};

#endif
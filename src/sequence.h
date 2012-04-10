//
//  sequence.h
//  express
//
//  Created by Adam Roberts on 1/24/12.
//  Copyright 2012 Adam Roberts. All rights reserved.
//

#ifndef express_sequence_h
#define express_sequence_h

#include <string>
#include "frequencymatrix.h"

/**
 * function to encode a nucleotide character to a size_t value
 * @param c the nucleotide character to be encoded
 * @return a size_t value encoding the nucleotide
 */
inline char ctoi(const char c)
{
    switch(c)
    {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return 0;
    }
}

inline char complement(const char c)
{
    switch(c)
    {
        case 0:
            return 3;
        case 1:
            return 2;
        case 2:
            return 1;
        case 3:
            return 0;
        default:
            assert(false);
            return 4;
    }    
}
/**
 * The Sequence class is used to store and access encoded nucleotide sequences.
 *  @author    Adam Roberts
 *  @date      2012
 *  @copyright Artistic License 2.0
 **/
class Sequence
{
    
public:
        
    /**
     * a member function that returns the encoded character at the given index
     * @param index the index of the encoded character to return (assumed to be < _len)
     * @return the encoded character at the given index
     */
    virtual size_t operator[](const size_t index) const = 0;
    
    virtual void update_obs(const size_t index, const size_t nuc, float mass) = 0;
    virtual void update_exp(const size_t index, const size_t nuc, float mass) = 0;
    virtual float get_prob(const size_t index, const size_t nuc) const = 0;
    virtual bool prob() const = 0;
    
    /**
     * a member function that returns the length of the encoded sequence
     * @return the length of the encoded sequence
     */
    virtual size_t length() const = 0;
    virtual bool empty() const = 0;
};

/**
 * The Sequence class is used to store and access encoded nucleotide sequences.
 *  @author    Adam Roberts
 *  @date      2012
 *  @copyright Artistic License 2.0
 **/
class SequenceFwd: public Sequence
{
    /**
     * a char array that stores the encoded sequence (not null terminated)
     */
    const char* _encoded_seq;
    
    FrequencyMatrix<float> _obs_seq;
    FrequencyMatrix<float> _exp_seq;
    
    size_t _len;
    
    bool _prob;
    
public:
    
    /**
     * dummy constructor
     */    
    SequenceFwd();
    
    /**
     * Sequence constructor encodes and stores the given nucleotide sequence
     * @param seq the nucleotide sequence to encode and store
     * @param rev a boolean if the sequence should be reverse complemented before encoding
     */   
    SequenceFwd(const std::string& seq, bool rev, bool prob=false);
    
    /**
     * Sequence copy constructor
     * @param other the Sequence object to copy
     */   
    SequenceFwd(const SequenceFwd& other);
    
    /**
     * Sequence assignment constructor
     * @param other the Sequence object to copy
     */   
    SequenceFwd& operator=(const SequenceFwd& other);
    
    /**
     * Sequence deconstructor. Deletes the char array.
     */   
    ~SequenceFwd();
    
    /**
     * a member function that encodes the given sequence and overwrites the current stored 
     * sequence with it
     * @param seq the nucleotide sequence to encode and store
     * @param rev a boolean if the sequence should be reverse complemented before encoding
     */   
    void set(const std::string& seq, bool rev);
    
    /**
     * a member function that returns the encoded character at the given index
     * @param index the index of the encoded character to return (assumed to be < _len)
     * @return the encoded character at the given index
     */
    size_t operator[](const size_t index) const;

    void update_obs(const size_t index, const size_t nuc, float mass);
    void update_exp(const size_t index, const size_t nuc, float mass);
    float get_prob(const size_t index, const size_t nuc) const;
    bool prob() const { return _prob; }
    bool empty() const { return _len==0; }
    size_t length() const { return _len; }
};

class SequenceRev: public Sequence
{

    SequenceFwd* _seq;
    
public:
    SequenceRev() {}
    SequenceRev(SequenceFwd& seq) : _seq(&seq) {}
    size_t length() const { return _seq->length(); }
    bool empty() const { return _seq->empty(); }
    void set(const std::string& seq, bool rev) { assert(false); }
    size_t operator[](const size_t index) const { return complement(_seq->operator[](length()-index-1)); }
    void update_obs(const size_t index, const size_t nuc, float mass) { _seq->update_obs(length()-index-1, complement(nuc), mass); }
    void update_exp(const size_t index, const size_t nuc, float mass) { _seq->update_exp(length()-index-1, complement(nuc), mass); }
    float get_prob(const size_t index, const size_t nuc) const { return _seq->get_prob(length()-index-1, complement(nuc)); }
    bool prob() const { return _seq->prob(); }
};

#endif

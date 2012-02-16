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

/**
 * The Sequence class is used to store and access encoded nucleotide sequences.
 *  @author    Adam Roberts
 *  @date      2012
 *  @copyright Artistic License 2.0
 **/
class Sequence
{
    /**
     * a char array that stores the encoded sequence (not null terminated)
     */
    const char* _encoded_seq;
    
    /**
     * a private size_t that stores the length of the sequence
     */
    size_t _len;
    
public:
    
    /**
     * dummy constructor
     */    
    Sequence();
    
    /**
     * Sequence constructor encodes and stores the given nucleotide sequence
     * @param seq the nucleotide sequence to encode and store
     * @param rev a boolean if the sequence should be reverse complemented before encoding
     */   
    Sequence(const std::string& seq, bool rev);
    
    /**
     * Sequence copy constructor
     * @param other the Sequence object to copy
     */   
    Sequence(const Sequence& other);
    
    /**
     * Sequence assignment constructor
     * @param other the Sequence object to copy
     */   
    Sequence& operator=(const Sequence& other);
    
    /**
     * Sequence deconstructor. Deletes the char array.
     */   
    ~Sequence();
    
    
    /**
     * a member function that encodes the given sequence and overwrites the current stored 
     * sequence with it
     * @param seq the nucleotide sequence to encode and store
     * @param rev a boolean if the sequence should be reverse complemented before encoding
     */   
    void set(const std::string& seq, bool rev);
    
    /**
     * a member function that returns true iff the encoded sequence has zero length
     * @return true iff the encoded sequence has zero length
     */
    bool empty() const { return _len==0; }
    

    /**
     * a member function that returns the encoded character at the given index
     * @param index the index of the encoded character to return (assumed to be < _len)
     * @return the encoded character at the given index
     */
    size_t operator[](const size_t index) const;

    /**
     * a member function that returns the length of the encoded sequence
     * @return the length of the encoded sequence
     */
    size_t length() const { return _len; }
};

#endif

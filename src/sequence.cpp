//
//  sequence.cpp
//  express
//
//  Created by Adam Roberts on 1/24/12.
//  Copyright 2012 UC Berkeley. All rights reserved.
//

#include "sequence.h"
#include <cassert>

using namespace std;

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

Sequence::Sequence():  _encoded_seq(NULL), _len(0) {}

Sequence::Sequence(const std::string& seq, bool rev) : _encoded_seq(NULL), _len(seq.length())
{
    set(seq, rev);
}

Sequence::Sequence(const Sequence& other) : _encoded_seq(NULL), _len(other.length())
{
    if (other._encoded_seq)
    {
        char* encoded_seq = new char[_len];
        std::copy(other._encoded_seq, other._encoded_seq + _len, encoded_seq);
        _encoded_seq = encoded_seq;
    }
}

Sequence& Sequence::operator=(const Sequence& other)
{
    if (other._encoded_seq)
    {
        _len = other.length();
        char* encoded_seq = new char[_len];
        std::copy(other._encoded_seq, other._encoded_seq + _len, encoded_seq);
        _encoded_seq = encoded_seq;
    }
    return *this;
}


Sequence::~Sequence()
{
    if (_encoded_seq)
        delete _encoded_seq;
}


void Sequence::set(const std::string& seq, bool rev)
{
    if (_encoded_seq)
        delete _encoded_seq;
    
    char* encoded_seq = new char[seq.length()]; 
    if (!rev)
    {
        for(size_t i = 0; i < seq.length(); i++)
        {
            encoded_seq[i] = ctoi(seq[i]);
        }
    }
    else
    {
        for(size_t i = 0; i < seq.length(); i++)
        {
            encoded_seq[seq.length()-1-i] = complement(ctoi(seq[i]));
        }
    }
    _encoded_seq = encoded_seq;
    _len = seq.length();
}

size_t Sequence::operator[](const size_t index) const
{
    assert(index < _len); 
    return _encoded_seq[index]; 
}

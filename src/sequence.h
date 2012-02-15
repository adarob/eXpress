//
//  sequence.h
//  express
//
//  Created by Adam Roberts on 1/24/12.
//  Copyright 2012 UC Berkeley. All rights reserved.
//

#ifndef express_sequence_h
#define express_sequence_h

#include <string>

//DOC
class Sequence
{
    const char* _encoded_seq;
    size_t _len;
    
public:
    
    Sequence();
    Sequence(const std::string& seq, bool rev);
    Sequence(const Sequence& other);
    Sequence& operator=(const Sequence& other);
    
    ~Sequence();
    
    void set(const std::string& seq, bool rev);
    bool empty() const { return _len==0; }
    size_t operator[](const size_t index) const;
    size_t length() const { return _len; }
};

#endif

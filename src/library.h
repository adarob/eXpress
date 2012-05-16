//
//  library.h
//  express
//
//  Created by Adam Roberts on 5/12/12.
//  Copyright (c) 2012 Adam Roberts. All rights reserved.
//

#ifndef express_library_h
#define express_library_h

#include <vector>

//DOC
/**
 * a struct for holding pointers to the global parameter tables (bias_table, mismatch_table, fld)
 */
struct Library
{
    MapParser* map_parser;
    FLD* fld;
    MismatchTable* mismatch_table;
    BiasBoss* bias_table;
    TargetTable* targ_table;
    size_t n;
    double mass_n;
    
    Library() : n(1), mass_n(0) {};
};

class Librarian
{
    std::vector<Library> _libs;
    size_t _curr;
    
public:
    Librarian(size_t num_libs): _libs(num_libs), _curr(0) {}
    Library& operator[](size_t i)
    {
        assert(i < _libs.size());
        return _libs[i];
    }
    
    const Library& curr_lib() const { return _libs[_curr]; }
    
    void set_curr(size_t i)
    {
        assert (i < _libs.size());
        _curr = i;
    }
    
    size_t size() const { return _libs.size(); }
};

#endif

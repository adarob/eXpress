//
//  fragments.h
//  expressionline
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 UC Berkeley. All rights reserved.
//

#ifndef FRAGMENTS_H
#define FRAGMENTS_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/thread.hpp>


struct FragMap
{
    std::string name;
    std::string ref;
    std::string seq_l;
    std::string seq_r;
    int left;
    int right;
    int mate_l;
    int length() const
    { 
        if (left >= 0 && right > 0)
            return right - left;
        return -1;
    }
    bool has_mate() { return mate_l >= 0; }
};

class Fragment
{
    std::vector<FragMap*> _frag_maps;
//    std::vector<FragMap> _open_mates;
    FragMap* _om;
    std::string _name;
    void add_open_mate(FragMap* om);
public:
    Fragment() { _om = NULL; }
    ~Fragment();
    bool add_map_end(FragMap* f);
    const std::string name() const { return _name; }
    const size_t num_maps() const { return _frag_maps.size(); }
    const std::vector<FragMap*>& maps() const { return _frag_maps; }
};

struct ParseThreadSafety
{
    Fragment* next_frag;
    boost::mutex proc_lk;
    boost::mutex parse_lk;
};

class MapParser
{
    std::ifstream _in;
    char* _line_buff;
    FragMap* _frag_buff;
    
    bool map_end_from_line();
public:
    MapParser(std::string sam_file);
    ~MapParser(){ delete _line_buff; }
    void spawn_thread(ParseThreadSafety* thread_safety);
    void threaded_parse(ParseThreadSafety* thread_safety);
    bool next_fragment(Fragment& f);
};

#endif
//
//  fragments.cpp
//  expressionline
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 UC Berkeley. All rights reserved.
//

#include "fragments.h"
#include <string.h>
#include <cstring>
using namespace std;


int cigar_length(const char* cigar_str)
{
    const char* p_cig = cigar_str;
    int tot_len = 0;
    while (*p_cig) 
    {
        char* t;
        int op_len = (int)strtol(p_cig, &t, 10);
        char op_char = toupper(*t);
        if (op_char == 'M' || op_char == 'N') 
        {
            tot_len += op_len;
        }
        
        p_cig = t + 1;
    }
    return tot_len;
}

char complement(char c)
{
    switch(c)
    {
        case 'a':
        case 'A':
            return 'T';
        case 'c':
        case 'C':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
            return 'A';
        default:
            return 'N';
    }
}

size_t reverse(size_t offset, size_t read_len)
{
    return read_len - offset - 1;
}


Fragment::~Fragment()
{
    for (size_t i = 0; i < num_maps(); i++)
    {
        delete _frag_maps[i];
    }
}

bool Fragment::add_map_end(FragMap* m)
{
    if (_name == "")
    {
        _name = m->name;
    }
    else if (_name != m->name)
    {
        return false;
    }
    
    add_open_mate(m);
    
    return true;
}

void Fragment::add_open_mate(FragMap* pm)
{
    if (!_om)
    {
        _om = pm;
    }
    else if (pm->ref == _om->ref && pm->mate_l == _om->left && _om->mate_l == pm->left)
    {
        if (pm->left < _om->left)
        {
            pm->right = _om->right;
            pm->seq_r = _om->seq_l;
            _frag_maps.push_back(pm);
            _om = NULL;
            delete _om;
        }
        else
        {
            _om->right = pm->right;
            _om->seq_r = pm->seq_l;
            _frag_maps.push_back(_om);
            delete pm;
        }
        _om = NULL;
    }
    else
    {
        delete _om;
        _om = pm;
    }
}
 
//void Fragment::add_open_mate(FragMap& om)
//{
//    for( vector<FragMap>::iterator it = _open_mates.begin(); it != _open_mates.end(); ++it)
//    {
//        FragMap& pm = *it;
//        if (pm.ref == om.ref && pm.mate_l == om.left && om.mate_l == pm.left)
//        {
//            if (pm.left < om.left)
//            {
//                pm.right = om.right;
//                _frag_maps.push_back(pm);
//            }
//            else
//            {
//                om.right = pm.right;
//                _frag_maps.push_back(om);
//            }
//            
//            found = true;
//            _open_mates.erase(it);
//            break;
//        }
//    }
//    
//    if (!found)
//        _open_mates.push_back(om);
//}

const size_t BUFF_SIZE = 9999;
const char* MISMATCH_TAG = "MD:Z:";

MapParser::MapParser(string sam_file)
{
    _line_buff = new char[BUFF_SIZE];
    _in.open(sam_file.c_str());
    if(!_in.is_open())
    {
        cerr << "Unable to open SAM file '" << sam_file << "'.\n" ; 
        exit(1);
    }
    _frag_buff = new FragMap();
    while(_in.good())
    {
        _in.getline(_line_buff, BUFF_SIZE-1, '\n');
        if (_line_buff[0] != '@')
        {
            map_end_from_line();
            break;
        }
    }
    
}

void MapParser::threaded_parse(ParseThreadSafety* thread_safety)
{
    ParseThreadSafety& ts = *thread_safety;
    bool fragments_remain = true;
    while (fragments_remain)
    {
        Fragment* frag = new Fragment(); 
        fragments_remain = next_fragment(*frag);
        ts.next_frag = frag;
        ts.proc_lk.unlock();
        ts.parse_lk.lock();
    }
    ts.next_frag = NULL;
    ts.proc_lk.unlock();
}

bool MapParser::next_fragment(Fragment& nf)
{
    nf.add_map_end(_frag_buff);    
    
    _frag_buff = new FragMap();
    
    while(_in.good())
    {        
        _in.getline (_line_buff, BUFF_SIZE-1, '\n');
        if (!map_end_from_line())
            continue;
        if (!nf.add_map_end(_frag_buff))
            break;
        _frag_buff = new FragMap();
    }
    
    return _in.good();
}

bool MapParser::map_end_from_line()
{
    FragMap& f = *_frag_buff;
    char *p = strtok(_line_buff, "\t");
    string read_seq;

    
    int i = 0;
    while (p && i <= 9) 
    {
        switch(i++)
        {
            case 0:
                f.name = p;
                break;
            case 1:
                break;
            case 2:
                if(p[0] == '*')
                    goto stop;
                f.ref = p;
                break;
            case 3:
                f.left = atoi(p)-1;
                break;
            case 4:
                break;
            case 5:
                f.right = f.left + cigar_length(p);
                break;
            case 6:
                if(p[0] != '=')
                    goto stop;
                break;
            case 7:
                f.mate_l = atoi(p)-1;
                break;
            case 8:
                break;
            case 9:
                f.seq_l = p;
                goto stop;
        }
        
        p = strtok(NULL, "\t");
    }
stop:
    return i >= 9;
}

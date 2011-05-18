//
//  fragments.cpp
//  expressionline
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 UC Berkeley. All rights reserved.
//

#include "fragments.h"
#include "transcripts.h"
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
            return 'T';
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
    
    for (size_t i = 0; i < _open_mates.size(); i++)
    {
        delete _open_mates[i];
    }
}

bool Fragment::add_map_end(FragMap* m)
{
    if (_name.empty())
    {
        _name = m->name;
    }
    else if (_name != m->name)
    {
        return false;
    }
    
    if (m->mate_l >= 0)
    {
        add_open_mate(m);
    }
    else
    {
        _frag_maps.push_back(m);
    }
    
    return true;
}

void Fragment::add_open_mate(FragMap* new_p)
{
    bool found = false;
    
    FragMap& nm = *new_p;
    for( vector<FragMap*>::iterator it = _open_mates.begin(); it != _open_mates.end(); ++it)
    {
        FragMap& om = **it;
        if (nm.trans_id == om.trans_id && nm.mate_l == om.left && om.mate_l == nm.left)
        {
            if (nm.left < om.left)
            {
                assert(nm.seq_r.empty());
                nm.right = om.right;
                nm.seq_r = om.seq_r;
            }
            else
            {
                assert(nm.seq_l.empty());
                nm.left = om.left;
                nm.seq_l = om.seq_l;
            }
            
            assert(nm.left_first == om.left_first);
            found = true;
            delete *it;
            _frag_maps.push_back(&nm);
            _open_mates.erase(it);
            break;
        }
    }
    
    if (!found)
        _open_mates.push_back(&nm);
}

//void Fragment::add_open_mate(FragMap* om)
//    {
//    if (!_om)
//    {
//        _om = pm;
//    }
//    else if (pm->trans_id == _om->trans_id && pm->mate_l == _om->left && _om->mate_l == pm->left)
//    {
//        if (pm->left < _om->left)
//        {
//            pm->right = _om->right;
//            pm->seq_r = _om->seq_r;
//            _frag_maps.push_back(pm);
//            _om = NULL;
//            delete _om;
//        }
//        else
//        {
//            _om->right = pm->right;
//            _om->seq_r = pm->seq_r;
//            _frag_maps.push_back(_om);
//            delete pm;
//        }
//        assert(_frag_maps.back().pair_status() == PAIRED);
//        _om = NULL;
//    }
//    else
//    {
//        delete _om;
//        _om = pm;
//    }
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
    int sam_flag;    
    bool paired;
    
    int i = 0;
    while (p && i <= 9) 
    {
        switch(i++)
        {
            case 0:
                f.name = p;
                break;
            case 1:
                sam_flag = atoi(p);
                paired = sam_flag & 0x1;
                f.left_first = (sam_flag & 0x40) && (sam_flag & 0x10) || !(sam_flag & 0x40) && !(sam_flag & 0x10);
                break;
            case 2:
                if(p[0] == '*')
                    goto stop;
                f.trans_id = hash_trans_name(p);
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
                if(paired && p[0] != '=')
                    goto stop;
                break;
            case 7:
                f.mate_l = atoi(p)-1;
                break;
            case 8:
                break;
            case 9:
                if (sam_flag & 0x10)
                    f.seq_r = p;
                else
                    f.seq_l = p;
                goto stop;
        }
        
        p = strtok(NULL, "\t");
    }
stop:
    return i == 9;
}

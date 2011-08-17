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

const size_t BUFF_SIZE = 9999;
const char* MISMATCH_TAG = "MD:Z:";

MapParser::MapParser(string sam_file)
{
    _line_buff = new char[BUFF_SIZE];
    if (sam_file.size() == 0)
    {
        _in = &cin;
    }
    else
    {
        ifstream* ifs = new ifstream(sam_file.c_str());
        if(!ifs->is_open())
        {
            cerr << "Unable to open SAM file '" << sam_file << "'.\n" ; 
            exit(1);
        } 
        _in = ifs;
    }

    _frag_buff = new FragMap();
    while(_in->good())
    {
        _in->getline(_line_buff, BUFF_SIZE-1, '\n');
        if (_line_buff[0] != '@')
        {
            map_end_from_line();
            break;
        }
    }
    
}

void MapParser::threaded_parse(ParseThreadSafety* thread_safety, TranscriptTable* trans_table)
{
    ParseThreadSafety& ts = *thread_safety;
    bool fragments_remain = true;
    while (fragments_remain)
    {
        Fragment * frag;
        while (fragments_remain)
        {
            frag = new Fragment(); 
            fragments_remain = next_fragment(*frag);
            if (frag->num_maps())
                break;
            delete frag;
            frag = NULL;
        }
        for (size_t i = 0; frag && i < frag->maps().size(); ++i)
        {
            FragMap& m = *(frag->maps()[i]);
            Transcript* t = trans_table->get_trans(m.trans_id);
            m.mapped_trans = t;
            assert(t->id() == m.trans_id);
        }
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
    
    while(_in->good())
    {        
        _in->getline (_line_buff, BUFF_SIZE-1, '\n');
        if (!map_end_from_line())
            continue;
        if (!nf.add_map_end(_frag_buff))
            break;
        _frag_buff = new FragMap();
    }
    
    return _in->good();
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
                if (sam_flag & 0x4)
                    goto stop;
                paired = sam_flag & 0x1;
                f.left_first = ((sam_flag & 0x40) && (sam_flag & 0x10)) || (!(sam_flag & 0x40) && !(sam_flag & 0x10));
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
    return i == 10;
}

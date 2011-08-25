//
//  mapparser.cpp
//  express
//
//  Created by Adam Roberts on 8/19/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "mapparser.h"
#include "main.h"
#include "fragments.h"
#include "transcripts.h"

using namespace std;

const size_t BUFF_SIZE = 9999;

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

int cigar_length(vector<BamTools::CigarOp> cigar_vec)
{
    int tot_len = 0;
    for (size_t i = 0; i < cigar_vec.size(); ++i)
    {
        if (cigar_vec[i].Type == 'M' || cigar_vec[i].Type == 'N') 
        {
            tot_len += cigar_vec[i].Length;
        }
    }
    
    return tot_len;
}

ThreadedMapParser::ThreadedMapParser(string file_name)
{
    if (file_name.size() == 0)
    {
        cout << "No alignment file specified. Expecting streaming input on stdin...\n";
        _parser = new SAMParser(&cin);
    }
    else
    {
        BamTools::BamReader* reader = new BamTools::BamReader();
        if (reader->Open(file_name))
        {
            _parser = new BAMParser(reader);
        }
        else
        {
            delete reader;
            cout << "Input is not in BAM format. Trying SAM...\n";
            ifstream* ifs = new ifstream(file_name.c_str());
            _parser = new SAMParser(ifs);
            if(!ifs->is_open())
            {
                cerr << "Unable to open SAM file '" << file_name << "'.\n" ; 
                exit(1);
            } 
        }
    }
}

void ThreadedMapParser::threaded_parse(ParseThreadSafety* thread_safety, TranscriptTable* trans_table)
{
    ParseThreadSafety& ts = *thread_safety;
    bool fragments_remain = true;
    while (fragments_remain)
    {
        Fragment * frag;
        while (fragments_remain)
        {
            frag = new Fragment(); 
            fragments_remain = _parser->next_fragment(*frag);
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

BAMParser::BAMParser(BamTools::BamReader* reader)
{
    _reader = reader;
    BamTools::BamAlignment a;
    
    // Get first valid FragMap
    _frag_buff = new FragMap();
    do {
        if (!_reader->GetNextAlignment(a))
        {
            cerr << "Error: Input BAM file contains no valid alignments." << endl;
        }
    } while(!map_end_from_alignment(a));
}

bool BAMParser::next_fragment(Fragment& nf)
{    
    nf.add_map_end(_frag_buff);
    
    _frag_buff = new FragMap();
    BamTools::BamAlignment a;
    
    while(true)
    {   
        if (!_reader->GetNextAlignment(a))
        {
            return false;
        }
        else if (!map_end_from_alignment(a))
        {
            continue;
        }
        else if (!nf.add_map_end(_frag_buff))
        {
            return true;
        }
        else
        {
            if (!_reader->GetNextAlignment(a))
                return false;
            _frag_buff = new FragMap();
        }
        
    }
}

bool BAMParser::map_end_from_alignment(BamTools::BamAlignment& a)
{
    FragMap& f = *_frag_buff;
    
    if (!a.IsMapped())
        return false;
    
    if (a.IsPaired() && (!a.IsMateMapped() || a.RefID != a.MateRefID))
        return false;
    
    f.left_first = (a.IsFirstMate() && !a.IsReverseStrand()) || (a.IsSecondMate() && a.IsReverseStrand());
    f.name = a.Name;
    f.trans_id = hash_trans_name(_reader->GetReferenceData()[a.RefID].RefName);
    f.left = a.Position;
    f.right = f.left + cigar_length(a.CigarData);
    f.mate_l = a.MatePosition;
    
    if (a.IsReverseStrand())
        f.seq_r = a.QueryBases;
    else
        f.seq_l = a.QueryBases;
    
    return true;
}

SAMParser::SAMParser(istream* in)
{
    _in = in;
    
    char line_buff[BUFF_SIZE];
    
    _frag_buff = new FragMap();
    while(_in->good())
    {
        _in->getline(line_buff, BUFF_SIZE-1, '\n');
        if (line_buff[0] != '@')
        {
            map_end_from_line(line_buff);
            break;
        }
    }
    
}

bool SAMParser::next_fragment(Fragment& nf)
{
    nf.add_map_end(_frag_buff);    
    
    _frag_buff = new FragMap();
    char line_buff[BUFF_SIZE];
    
    while(_in->good())
    {        
        _in->getline (line_buff, BUFF_SIZE-1, '\n');
        if (!map_end_from_line(line_buff))
            continue;
        if (!nf.add_map_end(_frag_buff))
            break;
        _frag_buff = new FragMap();
    }
    
    return _in->good();
}

bool SAMParser::map_end_from_line(char* line)
{
    FragMap& f = *_frag_buff;
    char *p = strtok(line, "\t");
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
                f.left_first = ((sam_flag & 0x40) && !(sam_flag & 0x10)) || (!(sam_flag & 0x40) && (sam_flag & 0x10));
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
                // skip if only one end mapped of paired-end read
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

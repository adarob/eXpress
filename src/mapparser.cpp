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

ThreadedMapParser::ThreadedMapParser(string in_file, string out_file) 
{
    bool is_sam = false;
    _writer = NULL;
    if (in_file.size() == 0)
    {
        cout << "No alignment file specified. Expecting streaming input on stdin...\n\n";
        _parser = new SAMParser(&cin);
        is_sam = true;
    }
    else
    {
        BamTools::BamReader* reader = new BamTools::BamReader();
        if (reader->Open(in_file))
        {
            _parser = new BAMParser(reader);
            if (out_file.size())
            {
                BamTools::BamWriter* writer = new BamTools::BamWriter();
                if (writer->Open(out_file, reader->GetHeader(), reader->GetReferenceData()))
                {
                    _writer = new BAMWriter(writer);
                }
                else
                {
                    cerr << "Unable to open output BAM file '" << out_file << "'.\n" ; 
                    exit(1);  
                }
            }
        }
        else
        {
            delete reader;
            cout << "Input is not in BAM format. Trying SAM...\n";
            ifstream* ifs = new ifstream(in_file.c_str());
            if(!ifs->is_open())
            {
                cerr << "Unable to open input SAM file '" << in_file << "'.\n" ; 
                exit(1);
            }
            _parser = new SAMParser(ifs);
            is_sam = true;
        }
    }
    
    if (is_sam && out_file.size())
    {
        ofstream* ofs = new ofstream(out_file.c_str());
        if(!ofs->is_open())
        {
            cerr << "Unable to open output SAM file '" << out_file << "'.\n" ; 
            exit(1);
        }
        *ofs << _parser->header();
        _writer = new SAMWriter(ofs);
    }
}

ThreadedMapParser::~ThreadedMapParser()
{
    delete _parser;
    if (_writer)
        delete _writer;
}

void ThreadedMapParser::threaded_parse(ParseThreadSafety* thread_safety, TranscriptTable* trans_table)
{
    ParseThreadSafety& ts = *thread_safety;
    bool fragments_remain = true;
    Fragment * last_frag = NULL;
    while (true)
    {
        Fragment* frag = NULL;
        while (running && fragments_remain)
        {
            frag = new Fragment(); 
            fragments_remain = _parser->next_fragment(*frag);
            if (frag->num_hits())
                break;
            delete frag;
            frag = NULL;
        }
        for (size_t i = 0; frag && i < frag->hits().size(); ++i)
        {
            FragHit& m = *(frag->hits()[i]);
            Transcript* t = trans_table->get_trans(m.trans_id);
            if (!t)
            {
                cerr << "ERROR: Target sequence not found. Input same FASTA file used in alignment.\n";
                exit(1);
            }
            m.mapped_trans = t;
            assert(t->id() == m.trans_id);
        }

        last_frag = ts.next_frag;
        
        {
            boost::unique_lock<boost::mutex> lock(ts.mut);
            if (running)
            {
                ts.next_frag = frag;
                ts.frag_clean = true;
                
                ts.cond.notify_one();

                while(ts.frag_clean)
                {
                    ts.cond.wait(lock);
                }
            }
        }

        if (last_frag && _writer)
            _writer->write_fragment(*last_frag);
        if (last_frag)
            delete last_frag;
        if (!frag)
            break;
        if (!running)
        {
            delete frag;
            break;
        }
    }
}

BAMParser::BAMParser(BamTools::BamReader* reader)
{
    _reader = reader;
    BamTools::BamAlignment a;
    
    foreach(const BamTools::RefData& ref, _reader->GetReferenceData())
    {
        _trans_index[ref.RefName] = _trans_index.size();
    }
    
    // Get first valid FragHit
    _frag_buff = new FragHit();
    do {
        if (!_reader->GetNextAlignment(a))
        {
            cerr << "ERROR: Input BAM file contains no valid alignments.\n";
            exit(1);
        }
    } while(!map_end_from_alignment(a));
}

bool BAMParser::next_fragment(Fragment& nf)
{    
    nf.add_map_end(_frag_buff);
    
    BamTools::BamAlignment a;
    _frag_buff = new FragHit();

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
        _frag_buff = new FragHit();
    }
}

bool BAMParser::map_end_from_alignment(BamTools::BamAlignment& a)
{
    FragHit& f = *_frag_buff;
    
    if (!a.IsMapped())
        return false;
    
    if (a.IsPaired() && (!a.IsMateMapped() || a.RefID != a.MateRefID))
        return false;
    
    f.left_first = (!a.IsPaired() && !a.IsReverseStrand()) || (a.IsFirstMate() && !a.IsReverseStrand()) || (a.IsSecondMate() && a.IsReverseStrand());
    f.name = a.Name;
    f.trans_id = a.RefID;
    f.left = a.Position;
    f.right = f.left + cigar_length(a.CigarData);
    f.mate_l = a.MatePosition;
    
    if (a.IsReverseStrand())
    {
        f.seq_r = a.QueryBases;
        f.bam_r = a;
    }
    else
    {
        f.seq_l = a.QueryBases;
        f.bam_l = a;
    }
    
    return true;
}

BAMWriter::BAMWriter(BamTools::BamWriter* writer) : _writer(writer) {}
BAMWriter::~BAMWriter() 
{ 
    _writer->Close(); 
    delete _writer;
}

void BAMWriter::write_fragment(Fragment& f)
{
    foreach(FragHit* hit, f.hits())
    {
        hit->bam_l.AddTag("PH","f",(float)hit->probability);
        hit->bam_r.AddTag("PH","f",(float)hit->probability);
        _writer->SaveAlignment(hit->bam_l);
        _writer->SaveAlignment(hit->bam_r);
    }
}

SAMParser::SAMParser(istream* in)
{
    _in = in;
    
    char line_buff[BUFF_SIZE];
    _frag_buff = new FragHit();
    _header = "";
    
    while(_in->good())
    {
        _in->getline(line_buff, BUFF_SIZE-1, '\n');
        if (line_buff[0] != '@')
        {
            break;
        }
        _header += line_buff;
        _header += "\n";
        
        string str(line_buff);
        size_t idx = str.find("SN:");
        if (idx!=string::npos)
        {
            str = str.substr(str.find("SN:")+3);
            str = str.substr(0,str.find_first_of("\n\t "));
            if (_trans_index.find(str) == _trans_index.end())
            {
                _trans_index[str] = _trans_index.size();
            }
            else
            {
                cerr << "Warning: Target '" << str << "' appears twice in the SAM index.\n";
            }
        }
    }

    while(!map_end_from_line(line_buff))
    { 
        if (!_in->good())
        {
            cerr << "ERROR: Input SAM file contains no valid alignments.\n";
            exit(1);
        }
        _in->getline(line_buff, BUFF_SIZE-1, '\n');
    }
    
}

bool SAMParser::next_fragment(Fragment& nf)
{
    nf.add_map_end(_frag_buff);    
    
    _frag_buff = new FragHit();
    char line_buff[BUFF_SIZE];
    
    while(_in->good())
    {        
        _in->getline (line_buff, BUFF_SIZE-1, '\n');
        if (!map_end_from_line(line_buff))
            continue;
        if (!nf.add_map_end(_frag_buff))
            break;
        _frag_buff = new FragHit();
    }
    
    return _in->good();
}

bool SAMParser::map_end_from_line(char* line)
{
    FragHit& f = *_frag_buff;
    string sam_line(line);
    char *p = strtok(line, "\t");
    int sam_flag=0;    
    bool paired=0;
    
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
                if (paired && (sam_flag & 0x8))
                    goto stop;
                f.left_first = ((sam_flag & 0x40) && !(sam_flag & 0x10)) || (!(sam_flag & 0x40) && (sam_flag & 0x10));
                break;
            case 2:
                if(p[0] == '*')
                    goto stop;
                f.trans_id = _trans_index[p];
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
                {
                    f.seq_r = p;
                    f.sam_r = sam_line;
                }
                else
                {
                    f.seq_l = p;
                    f.sam_l = sam_line;
                }
                goto stop;
        }
        
        p = strtok(NULL, "\t");
    }
stop:
    return i == 10;
}

SAMWriter::SAMWriter(ostream* out) : _out(out) {}
SAMWriter::~SAMWriter() 
{ 
    _out->flush(); 
    delete _out;
}

void SAMWriter::write_fragment(Fragment& f)
{
    foreach(FragHit* hit, f.hits())
    {
        *_out << hit->sam_l << " PH:f:" << (float)hit->probability << endl;
        *_out << hit->sam_r << " PH:f:" << (float)hit->probability << endl;
    }
}

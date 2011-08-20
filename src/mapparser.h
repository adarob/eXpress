//
//  mapparser.h
//  express
//
//  Created by Adam Roberts on 8/19/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#ifndef express_mapparser_h
#define express_mapparser_h

#include <vector>
#include <string>
#include <boost/thread.hpp>
#include <api/BamReader.h>

class Fragment;
class FragMap;
class TranscriptTable;

class Parser
{
public:
    
    /**
     * a member function that loads all mappings of the next fragment
     * @param f the empty Fragment to add mappings to
     * @return true if more reads remain in the SAM/BAM file/stream, false otherwise
     */
    virtual bool next_fragment(Fragment& f)=0;
};


/**
 * BAMParser class.
 *
 */
class BAMParser : public Parser
{
    BamTools::BamReader* _reader;
    
    BamTools::BamAlignment _align_buff;
    
    /**
     * a private pointer to the current fragment mapping being parsed
     */
    FragMap* _frag_buff;
    
    /**
     * a private function to parse a single alignment and store the data in _frag_buff
     * @return true if the mapping is valid and false otherwise
     */
    bool map_end_from_alignment(BamTools::BamAlignment& alignment);
    
public:
    
    /**
     * BAMParser constructor opens the file
     */
    BAMParser(BamTools::BamReader* reader);
    
    ~BAMParser() { delete _reader; }
    
    /**
     * a member function that loads all mappings of the next fragment
     * @param f the empty Fragment to add mappings to
     * @return true if more reads remain in the BAM file, false otherwise
     */
    bool next_fragment(Fragment& f);
};

/**
 * SAMParser class.  This class produces Fragment objects by parsing an input in SAM format.
 * The input may come from a file or stdin.
 */
class SAMParser : public Parser
{
    /**
     * a private pointer to the input stream (either stdin or file) in SAM format
     */
    std::istream* _in;
    
    /**
     * a private pointer to the current fragment mapping being parsed
     */
    FragMap* _frag_buff;
    
    /**
     * a private function to parse a single line and store the data in _frag_buff
     * @return true if the mapping is valid and false otherwise
     */
    bool map_end_from_line(char* line);
    
public:
    /**
     * SAMParser constructor removes the header, and parses the first line
     */
    SAMParser(std::istream* in);
    
    /**
     * a member function that loads all mappings of the next fragment
     * @param f the empty Fragment to add mappings to
     * @return true if more reads remain in the SAM file, false otherwise
     */
    bool next_fragment(Fragment& f);
};

/**
 * ParseThreadSafety struct.  This struct stores objects to allow for parsing to safely occur
 * on a separate thread from processing.
 */
struct ParseThreadSafety
{
    /**
     * a pointer to the next Fragment to be processed by the main thread
     */
    Fragment* next_frag;
    
    /**
     * a mutex to lock the main (processing) thread when next_frag has not yet been updated
     */
    boost::mutex proc_lk;
    
    /**
     * a mutex to lock the parsing thread when the next_frag pointer should be not modified
     */
    boost::mutex parse_lk;
};

/**
 * ThreadedMapParser class.
 */
class ThreadedMapParser
{
    Parser* _parser;
    
public:
    /**
     * ThreadedMapParser constructor determines what format the input is in and initializes the correct parser.
     */
    ThreadedMapParser(std::string input_file);
    
    ~ThreadedMapParser() { delete _parser; }
    /**
     * a member function that drives the parse thread
     * when all valid mappings of a fragment have been parsed, its mapped transcripts are found and the information
     * is passed in a Fragment object to the processing thread through the ParseThreadSafety struct
     * @param thread_safety a pointer to the struct containing shared locks and data with the processing thread
     * @param trans_table a pointer to the table of Transcript objects to lookup the mapped transcripts
     */
    void threaded_parse(ParseThreadSafety* thread_safety, TranscriptTable* trans_table);
    
};



#endif

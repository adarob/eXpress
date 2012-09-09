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
#include <boost/unordered_map.hpp>
#include <api/BamReader.h>
#include <api/BamWriter.h>

class ParseThreadSafety;
class Fragment;
struct FragHit;
class TargetTable;
struct Library;

typedef boost::unordered_map<std::string, size_t> TransIndex;

/**
 * The Parser class is an abstract class that can be a SAMParser or BAMParser.
 * It fills Fragment objects by parsing an input file in SAM/BAM format.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class Parser
{
public:
    virtual ~Parser(){};

    /**
     * a member function that returns a string version of the header
     * @return string version of the header
     */
    virtual const std::string header() const=0;
    
    /**
     * a member function that returns the target-to-index map
     * @return the target-to-index map
     */
    virtual const TransIndex& targ_index() const=0;
    
    /**
     * a member function that returns the target-to-length map
     * @return the target-to-length map
     */
    virtual const TransIndex& targ_lengths() const=0;
    
    /**
     * a member function that loads all mappings of the next fragment
     * @param f the empty Fragment to add mappings to
     * @return true if more reads remain in the SAM/BAM file/stream, false otherwise
     */
    virtual bool next_fragment(Fragment& f)=0;

    /**
     * a member function that resets the parser and rewinds to the beginning of the input
     */
    virtual void reset() = 0;
};

/**
 * The Writer class is an abstract class than can be a SAMWriter or BAMWriter.
 * It writes Fragment objects back to file (in SAM/BAM format) with per-mapping probabilistic assignments.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class Writer
{
    /**
     * a private bool that specifies if a single alignment should be sampled (true) or all output with
     * their respective probabilities (false)
     */
    bool _sample;
    
public:
    virtual ~Writer(){};
    
    /**
     * a member function that writes all mappings of the fragment to the ouptut
     * file along with their probabilities in the "XP" field
     * @param f the processed Fragment to output
     */
    virtual void write_fragment(Fragment& f)=0;
};


/**
 * The BAMParser class fills Fragment objects by parsing an input file in BAM format.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class BAMParser : public Parser
{
    /**
     * a private pointer to the BamReader object which direclty parses the BAM file
     */
    BamTools::BamReader* _reader;
    
    /**
     * the private target-to-index map
     */
    TransIndex _targ_index; 
    
    /**
     * the private target-to-length map 
     */
    TransIndex _targ_lengths;
    
    /**
     * a private pointer to the current fragment mapping being parsed
     */
    FragHit* _frag_buff;
    
    /**
     * a private function to parse a single alignment and store the data in _frag_buff
     * @return true if the mapping is valid and false otherwise
     */
    bool map_end_from_alignment(BamTools::BamAlignment& alignment);
    
public:
    
    /**
     * BAMParser constructor sets the reader
     */
    BAMParser(BamTools::BamReader* reader);
    
    /**
     * BAMParser destructor deletes the reader
     */
    ~BAMParser() { delete _reader; }
    
    /**
     * a member function that returns a string version of the header
     * @return string version of the header
     */
    const std::string header() const { return _reader->GetHeaderText(); }
    
    /**
     * a member function that returns the target-to-index map
     * @return the target-to-index map
     */
    const TransIndex& targ_index() const { return _targ_index; }
    
    /**
     * a member function that returns the target-to-length map
     * @return the target-to-length map
     */
    const TransIndex& targ_lengths() const { return _targ_lengths; }
    
    /**
     * a member function that loads all mappings of the next fragment
     * @param f the empty Fragment to add mappings to
     * @return true if more reads remain in the BAM file, false otherwise
     */
    bool next_fragment(Fragment& f);
    
    /**
     * a member function that resets the parser and rewinds to the beginning of the BAM file
     */
    void reset();
};

/**
 * The BAMWriter class writes Fragment objects back to file (in BAM format) with per-mapping probabilistic assignments.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class BAMWriter : public Writer
{
    /**
     * a private pointer to the BamWriter object which direclty writes to the BAM file
     */
    BamTools::BamWriter* _writer;
    
    /**
     * a private bool that specifies if a single alignment should be sampled (true) or all output with
     * their respective probabilities (false)
     */
    bool _sample;
    
public:
    
    /**
     * BAMWriter constructor stores a pointer to the BAM file
     * @param writer pointer to the BAM file object
     * @param sample specifies if a single alignment should be sampled (true) or all output with
     * their respective probabilities (false)
     */
    BAMWriter(BamTools::BamWriter* writer, bool sample);
    
    /**
     * BAMWriter destructor flushes and deletes the Bamtools::BamWriter
     */
    ~BAMWriter();
    
    /**
     * a member function that writes all mappings of the fragment to the ouptut
     * file in BAM format along with their probabilities in the "XP" field
     * @param f the processed Fragment to output
     */
    void write_fragment(Fragment& f);
};


/**
 * The SAMParser class fills Fragment objects by parsing an input in SAM format.
 * The input may come from a file or stdin.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class SAMParser : public Parser
{
    /**
     * a private pointer to the input stream (either stdin or file) in SAM format
     */
    std::istream* _in;
    
    /**
     * the private target-to-index map
     */
    TransIndex _targ_index; 
    
    /**
     * the private target-to-length map
     */
    TransIndex _targ_lengths;
    
    /**
     * a private pointer to the current fragment mapping being parsed
     */
    FragHit* _frag_buff;
    
    /**
     * a private string storing the SAM header
     */
    std::string _header;
    
    /**
     * a private function to parse a single line and store the data in _frag_buff
     * @return true if the mapping is valid and false otherwise
     */
    bool map_end_from_line(char* line);
    
    /**
     * a member function that resets the parser and rewinds to the beginning of the SAM file
     */
    void reset();
    
public:
    /**
     * SAMParser constructor removes the header, and parses the first line
     */
    SAMParser(std::istream* in);
    
    /**
     * a member function that returns a string version of the header
     * @return string version of the header
     */
    const std::string header() const { return _header; }
    
    /**
     * a member function that returns the target-to-index map
     * @return the target-to-index map
     */
    const TransIndex& targ_index() const { return _targ_index; }
    
    /**
     * a member function that returns the target-to-length map
     * @return the target-to-length map
     */
    const TransIndex& targ_lengths() const { return _targ_lengths; }
    
    /**
     * a member function that loads all mappings of the next fragment
     * when the next fragment is reached, the current alignment is left in the
     * _frag_buff for the next call
     * @param f the empty Fragment to add mappings to
     * @return true if more reads remain in the SAM file, false otherwise
     */
    bool next_fragment(Fragment& f);
};

/**
 * The SAMWriter class writes Fragment objects back to file (in SAM format) with per-mapping probabilistic assignments.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class SAMWriter : public Writer
{
    /**
     * a private pointer to the output stream to which the alignments are written in SAM format
     */
    std::ostream* _out;
    
    /**
     * a private bool that specifies if a single alignment should be sampled (true) or all output with
     * their respective probabilities (false)
     */
    bool _sample;
    
public:
    
    /**
     * SAMWriter constructor stores a pointer to the output stream
     * @param out SAM output stream
     * @param sample specifies if a single alignment should be sampled (true) or all output with
     * their respective probabilities (false)
     */
    SAMWriter(std::ostream* out, bool sample);

    /**
     * SAMWriter destructor flushes and deletes output stream
     */
    ~SAMWriter();
    
    /**
     * a member function that writes all mappings of the fragment to the ouptut
     * file in SAM format along with their probabilities in the "XP" field
     * @param f the processed Fragment to output
     */
    void write_fragment(Fragment& f);
};


/**
 * The MapParser class is meant to be run on as a separate thread from the main processing.
 * Once started, this thread will read input from a file or stream in SAM/BAM format, parse, and collect 
 * read alignments into fragment alignments, and fragment alignments into fragments, which are placed on a
 * buffer for the processing thread.  Once the processing thread copies the fragment address from the buffer,
 * the parser is unlocked to load the next fragment.  The process stops when EOF is reached
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class MapParser
{
    /**
     * a private pointer to the Parser object that will read the input in SAM/BAM format
     */
    Parser* _parser;
    
    /**
     * a private pointer to the Writer object that will write the output in SAM/BAM format
     */
    Writer* _writer;
    
    //DOC
    std::string _in_file;
    
    Library* _lib;
    
    /**
     * a private boolean specifying whether to output the modified alignments
     */
    bool _write_active;
    
public:
    /**
     * MapParser constructor determines what format the input is in and initializes the correct parser.
     * @param input_file string containing the path to the input SAM/BAM file
     * @param output_file string containing the path the output file less its extension (empty if writing is to be disabled)
     */
    MapParser(std::string input_file, std::string output_file, Library* lib, bool write_active);
    
    /**
     * MapParser destructor deletes the parser and writer (if it exists).
     */
    ~MapParser();
    
    //DOC
    std::string& in_file_name() { return _in_file; }
    
    /**
     * a member function that drives the parse thread
     * when all valid mappings of a fragment have been parsed, its mapped targets are found and the information
     * is passed in a Fragment object to the processing thread through the ParseThreadSafety struct
     * @param thread_safety a pointer to the struct containing shared queues with the processing thread
     * @param stop_at a size_t indicating how many reads to process before stopping (disabled if 0)
     */
    void threaded_parse(ParseThreadSafety* thread_safety, size_t stop_at=0, size_t num_neighbors=0);
    
    /**
     * a member function that returns the target-to-index map
     * @return the target-to-index map
     */
    const TransIndex& targ_index() { return _parser->targ_index(); }
    
    /**
     * a member function that returns the target-to-length map
     * @return the target-to-length map
     */
    const TransIndex& targ_lengths() { return _parser->targ_lengths(); }
    
    /**
     * a member function that sets the write-active status of the parser
     * this specifies whether or not the alignments (sampled or with probs) sould be ouptut
     * @param b updated write-active status
     */    
    void write_active(bool b) { _write_active = b; }
    
    /**
     * a member function that resets the input parser
     */
    void reset_reader() { _parser->reset(); }
    
};



#endif

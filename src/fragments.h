//
//  fragments.h
//  express
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#ifndef FRAGMENTS_H
#define FRAGMENTS_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>

typedef size_t TransID;
class Transcript;
class TranscriptTable;

/**
 * PairStatus enum
 * PAIRED denotes that the fragment mapping has both paired end reads
 * LEFT_ONLY denotes that the single read is not reverse complemented => its left end is the left fragment end
 * RIGHT_ONLY denotes that the single read is reverse complemented => its right end is the right fragment end
 */
enum PairStatus { PAIRED, LEFT_ONLY, RIGHT_ONLY };

/**
 *  The FragMap struct stores the information for a single (multi-)mapping of a fragment.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
struct FragMap
{
    /**
     * a public string for the SAM "Query Template Name" (fragment name)
     */
    std::string name;
    
    /**
     * a public TransID for the transcript mapped to
     */
    TransID trans_id;
    
    /**
     * a public pointer to the transcript mapped to
     */
    Transcript* mapped_trans;
    
    /**
     * a public string containing the "left" read sequence (first according to SAM flag)
     */
    std::string seq_l;
    
    /**
     * a public string containing the "right" read sequence (second according to SAM flag)
     */
    std::string seq_r;
    
    /**
     * a public int containing the 0-based leftmost coordinate mapped to in the transcript
     * valid only if PairStatus is PAIRED or LEFT_ONLY
     */
    int left;
    
    /**
     * a public int containing the position following the 0-based rightmost coordinate mapped to in the transcript
     * valid only if PairStatus is PAIRED or RIGHT_ONLY
     */
    int right;
    
    /**
     * a public int containing the left position for the mate of the first read read in from the SAM file 
     * 0 if single-end fragment
     * this is temporarily used to help find the mate, but is not important later on
     */
    int mate_l;
    
    /**
     * a public bool specifying that the "right" (second according to SAM flag) is reverse complemented when true
     * and the "left" (first according to SAM flag) is reverse complemented when false
     * in other words, the "left" read is truly left of the "right" read in transcript coordinate space when true
     */
    bool left_first;
    
    /**
     * a member function returning the length of the fragment according to this mapping
     * note, that this result will be invalid if the fragment is single-end
     * @return int length of fragment mapping
     */
    int length() const
    {
        assert(pair_status() == PAIRED);
        return right - left;
    }
    
    /**
     * a member function returning whether the mapping is PAIRED, LEFT_ONLY, or RIGHT_ONLY
     * LEFT_ONLY denotes that the single read is not reverse complemented => its left end is the left fragment end
     * RIGHT_ONLY denotes that the single read is reverse complemented => its right end is the right fragment end
     * @return PairStatus the pair status of the mapping
     */
    PairStatus pair_status() const
    { 
        if (!seq_l.empty() && !seq_r.empty())
            return PAIRED;
        if (seq_l.empty())
            return LEFT_ONLY;
        return RIGHT_ONLY;
    }
};

/**
 * The Fragment class stores information for all multi-mappings of a single fragment.
 * By design, only paired-end mappings of paired-end reads will be accepted.  All mappings
 * of single-end reads will be accepted.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
class Fragment
{
    
    /**
     * a private vector of FragMap pointers containing all multi-mappings of the fragment
     */
    std::vector<FragMap*> _frag_maps;
    
    /**
     * a private vector of FragMap pointers containing single-end mappings whose pairs have not been found
     * this is temporarily used when parsing the file to find mates, but is not important later on
     */
    std::vector<FragMap*> _open_mates;
    
    /**
     * a private string for the SAM "Query Template Name" (fragment name)
     */
    std::string _name;
    
    /**
     * a private function that searches for the mate of the given read mapping
     * if found, the mates are combined into a single fragment and added to _frag_maps
     * if not found, the read mapping is added to open_mates
     */
    void add_open_mate(FragMap* om);

public:
    
    /**
     * Fragment destructor deletes all FragMap objects pointed to by the Fragment
     */
    ~Fragment();
    
    /**
     * a member function that adds a new FragMap (single read at this point) to the Fragment
     * if it is the first FragMap, it sets the Fragment name and is added to _open_mates,
     * if the fragment is not paired, it is added to _frag_maps,
     * otherwise, add_open_mate is called
     * @param f the FragMap to be added
     */
    bool add_map_end(FragMap* f);
    
    /**
     * a member function that returns the SAM "Query Template Name" (fragment name)
     * @return the string SAM "Query Template Name" (fragment name)
     */
    const std::string name() const { return _name; }
    
    /**
     * a member function that returns the number of multi-mappings for the fragment
     * @return number of multi-mappings for fragment
     */
    const size_t num_maps() const { return _frag_maps.size(); }
    
    /**
     * a member function that returns FragMap multi-mappings of the fragment
     * @return a vector containing pointers to the FragMap multi-mappings
     */
    const std::vector<FragMap*>& maps() const { return _frag_maps; }
};


#endif
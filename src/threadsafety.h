//
//  thread_safety.h
//  express
//
//  Created by Adam Roberts on 1/25/12.
//  Copyright 2012 Adam Roberts. All rights reserved.
//

#ifndef express_thread_safety_h
#define express_thread_safety_h

#include <boost/thread.hpp>

class Fragment;

/**
 * The ParseThreadSafety struct stores objects to allow for parsing to safely occur
 * on a separate thread from processing.
 *  @author    Adam Roberts
 *  @date      2011
 *  @copyright Artistic License 2.0
 **/
struct ParseThreadSafety
{
    /**
     * a pointer to the next Fragment to be processed by the main thread
     */
    Fragment* next_frag;
    
    /**
     * a mutex for the conditional variable
     */
    boost::mutex mut;
    
    /**
     * a conditional variable where the processor waits for a new Fragment and the parser waits
     * for the Fragment pointer to be copied by the processor
     */
    boost::condition_variable cond;
    
    /**
     * a bool specifying the condition that the current next_frag pointer is clean, meaning that it
     * hasn't been copied by the processor
     */
    bool frag_clean;
    
    ParseThreadSafety() : next_frag(NULL), frag_clean(false) {}
};

#endif

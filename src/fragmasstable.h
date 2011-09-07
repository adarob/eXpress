//
//  fragmasstable.h
//  express
//
//  Created by Adam Roberts on 9/7/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#ifndef express_fragmasstable_h
#define express_fragmasstable_h

#include <vector>

const size_t INT_BIT = sizeof(int)*CHAR_BIT; 

class FragMassTable
{
    double _ff_param;
    
    double _n;
    double _mass_n;
    double _cum_mass_n;
    
    std::vector<double> _n_table;
    std::vector<double> _mass_table;
    std::vector<double> _cum_mass_table;
    
    size_t n_to_index (size_t n);
public:
    FragMassTable(double ff_param);
    
    double next_frag_mass();
    double next_frag_mass(size_t n, double curr_mass);
    
    /**
     * 
     * @param n number of fragments observed
     * @param mass mass of next fragment for returned number of fragments
     * @param cum_mass cumulative mass for returned number of fragments
     * @return nearest (floor) number of fragments for stored values
     */
    size_t nearest_stored_mass(size_t n, double& mass, double& cum_mass);
};

#endif

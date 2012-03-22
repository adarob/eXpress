//
//  markovmodel.h
//  express
//
//  Created by Adam Roberts on 1/24/12.
//  Copyright 2012 Adam Roberts. All rights reserved.
//

#ifndef express_markovmodel_h
#define express_markovmodel_h

#include <vector>
#include <string>

class Sequence;
class FrequencyMatrix;

class MarkovModel
{
    int _order;
    int _window_size;
    int _num_pos;
    std::vector<FrequencyMatrix> _params;
    bool _logged;
    
public:
    MarkovModel(size_t order, size_t window_size, size_t num_pos, double alpha, bool logged=true);
    
    double transition_prob(size_t p, size_t cond, size_t curr) const;

    double seq_prob(const Sequence& seq, int left) const;
        
    double marginal_prob(size_t w, size_t nuc) const;

    void update(const Sequence& seq, int left, double mass);
    
    void fast_learn(const Sequence& seq, double mass, const std::vector<double>& fl_cdf);
    void calc_marginals();
    
    void set_logged(bool logged);
};


#endif

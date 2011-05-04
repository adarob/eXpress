//
//  main.cpp
//  expressionline2
//
//  Created by Adam Roberts on 3/23/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#define ENABLE_THREADS 0

#include <boost/thread.hpp>
#include <cassert>
#include <iostream>
#include <cmath>
#include "transcripts.h"
#include "fld.h"
#include "fragments.h"
#include "biascorrection.h"
#include "mismatchmodel.h"

using namespace std;

double expr_alpha = .001;
double fld_alpha = 1;
double bias_alpha = 1;
double mm_alpha = 1;
int max_len = 800;

boost::mutex proc_lock;

void unlock_proc()
{
    proc_lock.unlock();
}

double log_sum(double x, double y)
{
    return x+log(1+exp(y-x));
}

void process_fragment(Fragment* frag_p, TranscriptTable* trans_table, FLD* fld, BiasBoss* bias_table, MismatchTable* mismatch_table)
{
#if ENABLE_THREADS
   	boost::this_thread::at_thread_exit(unlock_proc);
#endif
    Fragment& frag = *frag_p;
    
    if (frag.num_maps()==0)
    {
        delete frag_p;
        return;
    }
    if (frag.num_maps()==1)
    {
        const FragMap& m = *frag.maps()[0];
        Transcript* t  = trans_table->get_trans(m.ref);
        if (!t)
        {
            cerr << "ERROR: Transcript " << m.ref << " not found in reference fasta.";
            exit(1);
        }
        double mass = 1.0;
        t->add_mass(mass);
        fld->add_val(m.length(), mass);
        if (bias_table)
            bias_table->update_observed(m, *t, mass);
        mismatch_table->update(m, *t, mass);
        delete frag_p;
        return;
    }
    
    vector<double> likelihoods(frag.num_maps());
    vector<Transcript*> transcripts(frag.num_maps());
    double total_likelihood = 0.0;
    for(size_t i = 0; i < frag.num_maps(); ++i)
    {
        const FragMap& m = *frag.maps()[i];
        Transcript* t = trans_table->get_trans(m.ref);
        transcripts[i] = t;
        assert(t->id() == m.ref);
        //t->update_transcript_bias();
        likelihoods[i] = t->log_likelihood(m);
        if (!likelihoods[i])
            continue;
        total_likelihood = (i) ? log_sum(total_likelihood, likelihoods[i]):likelihoods[i];
    }
    if (total_likelihood == 0)
    {
        delete frag_p;
        return;
    }
    double total_mass = 0.0;
    for(size_t i = 0; i < frag.num_maps(); ++i)
    {
        if (!likelihoods[i])
            continue;
        const FragMap& m = *frag.maps()[i];
        Transcript* t  = transcripts[i];
        double mass = exp(likelihoods[i]-total_likelihood);
        assert(!isnan(mass));
        t->add_mass(mass);
        if (bias_table)
            bias_table->update_observed(m, *t, mass);
        mismatch_table->update(m, *t, mass);
        fld->add_val(m.length(), mass);
        total_mass += mass;
    }
    assert(abs(total_mass-1.0) < .0001);
    delete frag_p;
}

void threaded_calc_abundances(MapParser& map_parser, TranscriptTable* trans_table, FLD* fld, BiasBoss* bias_table, MismatchTable* mismatch_table)
{
    ParseThreadSafety ts;
    ts.proc_lk.lock();
    ts.parse_lk.lock();
    boost::thread parse(&MapParser::threaded_parse, &map_parser, &ts);
    size_t count = 0;
    while(true)
    {
        ts.proc_lk.lock();
        Fragment* frag = ts.next_frag;
        ts.parse_lk.unlock();
        if (!frag)
        {
            ts.proc_lk.unlock();
            return;
        }
        count += frag->num_maps();
        if (count % 100000 < frag->num_maps() )
            cout<< count << "\n";
        process_fragment(frag, trans_table, fld, bias_table, mismatch_table);
    }
}

void calc_abundances(MapParser& map_parser, TranscriptTable* trans_table, FLD* fld, BiasBoss* bias_table, MismatchTable* mismatch_table)
{
    bool fragments_remain = true;
    size_t count = 0;
    while(fragments_remain)
    {
        Fragment* frag = new Fragment();
        fragments_remain = map_parser.next_fragment(*frag);
        count += frag->num_maps();
        if (count % 10000 < frag->num_maps() )
            cout<< count << "\n";
#if ENABLE_THREADS
        proc_lock.lock();
        boost::thread process(process_fragment, frag, trans_table, fld, bias_table, mismatch_table);
#else
        process_fragment(frag, trans_table, fld, bias_table, mismatch_table);
#endif
    }
}

int main (int argc, const char * argv[])
{
    FLD fld(fld_alpha, max_len);
    BiasBoss bias_table(bias_alpha);
    MismatchTable mismatch_table(mm_alpha);
    MapParser map_parser("/Users/adarob/projects/Tuxedo/devel/expressionline2/test_data/hits.sam");
    BiasBoss* b= &bias_table;
    //b = NULL;
    TranscriptTable trans_table("/Users/adarob/projects/Tuxedo/devel/expressionline2/test_data/mtg_out.fa", expr_alpha, &fld, b, &mismatch_table);
    trans_table.threaded_bias_update();
    boost::thread bias_update(&TranscriptTable::threaded_bias_update, &trans_table);
    threaded_calc_abundances(map_parser, &trans_table, &fld, b, &mismatch_table);
    mismatch_table.output("/Users/adarob/projects/Tuxedo/devel/expressionline2/test_data");
    trans_table.output_expression("/Users/adarob/projects/Tuxedo/devel/expressionline2/test_data");
    return 0;
}

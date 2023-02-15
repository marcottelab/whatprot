/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "bootstrap-fit.h"

// Standard C++ library headers:
#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "fitters/hmm-fitter.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::cout;
using std::default_random_engine;
using std::sort;
using std::uniform_int_distribution;
using std::vector;
}  // namespace

void bootstrap_fit(unsigned int num_timesteps,
                   unsigned int num_channels,
                   double stopping_threshold,
                   const SequencingModel& seq_model,
                   const SequencingSettings& seq_settings,
                   const DyeSeq& dye_seq,
                   const vector<Radiometry>& radiometries,
                   unsigned int num_bootstrap_rounds,
                   double confidence_interval,
                   vector<SequencingModel>* seq_models) {
    HMMFitter fitter(num_timesteps,
                     num_channels,
                     stopping_threshold,
                     seq_model,
                     seq_settings,
                     dye_seq);
    // This prng is not thread safe and having one copy for each thread risks
    // issues with randomness. May be worth exploring multi-threaded randomness
    // options.
    default_random_engine generator(time_based_seed());
    // The uniform_int_distribution uses both limits inclusively!
    uniform_int_distribution<> rand_idx(0, radiometries.size() - 1);
    seq_models->resize(num_bootstrap_rounds);
#pragma omp parallel for schedule(dynamic, 1)
    for (unsigned int i = 0; i < num_bootstrap_rounds; i++) {
        vector<const Radiometry*> subsample;
#pragma omp critical
        for (unsigned int j = 0; j < radiometries.size(); j++) {
            subsample.push_back(&radiometries[rand_idx(generator)]);
        }
        // end pragma omp critical
        (*seq_models)[i] = fitter.fit(subsample);
    }
    // end pragma omp parallel for
    unsigned int ci_min_idx = (unsigned int)((1 - confidence_interval) / 2
                                             * (double)num_bootstrap_rounds);
    unsigned int ci_max_idx =
            (unsigned int)(((1 - (1 - confidence_interval) / 2))
                           * (double)num_bootstrap_rounds);
    SequencingModel ci_min(num_channels);
    SequencingModel ci_max(num_channels);
    sort(seq_models->begin(),
         seq_models->end(),
         [](SequencingModel a, SequencingModel b) -> bool {
             return a.p_edman_failure < b.p_edman_failure;
         });
    ci_min.p_edman_failure = (*seq_models)[ci_min_idx].p_edman_failure;
    ci_max.p_edman_failure = (*seq_models)[ci_max_idx].p_edman_failure;
    sort(seq_models->begin(),
         seq_models->end(),
         [](SequencingModel a, SequencingModel b) -> bool {
             return a.p_initial_detach < b.p_initial_detach;
         });
    ci_min.p_initial_detach = (*seq_models)[ci_min_idx].p_initial_detach;
    ci_max.p_initial_detach = (*seq_models)[ci_max_idx].p_initial_detach;
    sort(seq_models->begin(),
         seq_models->end(),
         [](SequencingModel a, SequencingModel b) -> bool {
             return a.p_cyclic_detach < b.p_cyclic_detach;
         });
    ci_min.p_cyclic_detach = (*seq_models)[ci_min_idx].p_cyclic_detach;
    ci_max.p_cyclic_detach = (*seq_models)[ci_max_idx].p_cyclic_detach;
    sort(seq_models->begin(),
         seq_models->end(),
         [](SequencingModel a, SequencingModel b) -> bool {
             return a.p_initial_break_n < b.p_initial_break_n;
         });
    ci_min.p_initial_break_n = (*seq_models)[ci_min_idx].p_initial_break_n;
    ci_max.p_initial_break_n = (*seq_models)[ci_max_idx].p_initial_break_n;
    sort(seq_models->begin(),
         seq_models->end(),
         [](SequencingModel a, SequencingModel b) -> bool {
             return a.p_cyclic_break_n < b.p_cyclic_break_n;
         });
    ci_min.p_cyclic_break_n = (*seq_models)[ci_min_idx].p_cyclic_break_n;
    ci_max.p_cyclic_break_n = (*seq_models)[ci_max_idx].p_cyclic_break_n;
    for (unsigned int c = 0; c < num_channels; c++) {
        sort(seq_models->begin(),
             seq_models->end(),
             [c](SequencingModel a, SequencingModel b) -> bool {
                 return a.channel_models[c]->p_initial_bleach
                        < b.channel_models[c]->p_initial_bleach;
             });
        ci_min.channel_models[c]->p_initial_bleach =
                (*seq_models)[ci_min_idx].channel_models[c]->p_initial_bleach;
        ci_max.channel_models[c]->p_initial_bleach =
                (*seq_models)[ci_max_idx].channel_models[c]->p_initial_bleach;
        sort(seq_models->begin(),
             seq_models->end(),
             [c](SequencingModel a, SequencingModel b) -> bool {
                 return a.channel_models[c]->p_cyclic_bleach
                        < b.channel_models[c]->p_cyclic_bleach;
             });
        ci_min.channel_models[c]->p_cyclic_bleach =
                (*seq_models)[ci_min_idx].channel_models[c]->p_cyclic_bleach;
        ci_max.channel_models[c]->p_cyclic_bleach =
                (*seq_models)[ci_max_idx].channel_models[c]->p_cyclic_bleach;
        sort(seq_models->begin(),
             seq_models->end(),
             [c](SequencingModel a, SequencingModel b) -> bool {
                 return a.channel_models[c]->p_dud < b.channel_models[c]->p_dud;
             });
        ci_min.channel_models[c]->p_dud =
                (*seq_models)[ci_min_idx].channel_models[c]->p_dud;
        ci_max.channel_models[c]->p_dud =
                (*seq_models)[ci_max_idx].channel_models[c]->p_dud;
    }
    cout << "lower bounds: " << ci_min.debug_string() << "\n";
    cout << "upper bounds: " << ci_max.debug_string() << "\n";
}

}  // namespace whatprot

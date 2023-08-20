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

double bootstrap_fit(unsigned int num_timesteps,
                     unsigned int num_channels,
                     double stopping_threshold,
                     double max_runtime,
                     const SequencingModel& seq_model,
                     const SequencingSettings& seq_settings,
                     const DyeSeq& dye_seq,
                     const vector<Radiometry>& radiometries,
                     unsigned int num_bootstrap_rounds,
                     double confidence_interval,
                     vector<SequencingModel>* seq_models,
                     vector<double>* log_ls,
                     SequencingModel* best_seq_model,
                     double* step_size) {
    HMMFitter fitter(num_timesteps,
                     num_channels,
                     stopping_threshold,
                     max_runtime,
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
    log_ls->resize(num_bootstrap_rounds);
    double best_log_l = 0.0;
    double worst_bootstrap_step_size;
#pragma omp parallel for schedule(dynamic, 1)
    for (unsigned int i = 0; i <= num_bootstrap_rounds; i++) {
        // Extra case so we can get base results to display. Otherwise we
        // subsample.
        if (i == num_bootstrap_rounds) {
            best_log_l = fitter.fit(radiometries, best_seq_model, step_size);
        } else {
            vector<const Radiometry*> subsample;
#pragma omp critical
            for (unsigned int j = 0; j < radiometries.size(); j++) {
                subsample.push_back(&radiometries[rand_idx(generator)]);
            }
            // end pragma omp critical
            double bootstrap_step_size;
            (*log_ls)[i] = fitter.fit(
                    subsample, &(*seq_models)[i], &bootstrap_step_size);
            // This should probably be a reduction instead of an omp critical
            // but I'm lazy and given the omp critical above this probably has a
            // marginal negative impact on runtime.
#pragma omp critical
            if (bootstrap_step_size > worst_bootstrap_step_size) {
                worst_bootstrap_step_size = bootstrap_step_size;
            }
        }
    }
    // end pragma omp parallel for
    // Temporary copy of seq_models so we can sort without messing up
    // concordance with log_ls.
    vector<SequencingModel> smz = *seq_models;
    unsigned int ci_min_idx = (unsigned int)((1 - confidence_interval) / 2
                                             * (double)num_bootstrap_rounds);
    unsigned int ci_max_idx =
            (unsigned int)(((1 - (1 - confidence_interval) / 2))
                           * (double)num_bootstrap_rounds);
    SequencingModel ci_min(num_channels);
    SequencingModel ci_max(num_channels);
    sort(smz.begin(),
         smz.end(),
         [](SequencingModel a, SequencingModel b) -> bool {
             return a.p_edman_failure < b.p_edman_failure;
         });
    ci_min.p_edman_failure = smz[ci_min_idx].p_edman_failure;
    ci_max.p_edman_failure = smz[ci_max_idx].p_edman_failure;
    sort(smz.begin(),
         smz.end(),
         [](SequencingModel a, SequencingModel b) -> bool {
             return a.p_detach < b.p_detach;
         });
    ci_min.p_detach = smz[ci_min_idx].p_detach;
    ci_max.p_detach = smz[ci_max_idx].p_detach;
    sort(smz.begin(),
         smz.end(),
         [](SequencingModel a, SequencingModel b) -> bool {
             return a.p_initial_block < b.p_initial_block;
         });
    ci_min.p_initial_block = smz[ci_min_idx].p_initial_block;
    ci_max.p_initial_block = smz[ci_max_idx].p_initial_block;
    sort(smz.begin(),
         smz.end(),
         [](SequencingModel a, SequencingModel b) -> bool {
             return a.p_cyclic_block < b.p_cyclic_block;
         });
    ci_min.p_cyclic_block = smz[ci_min_idx].p_cyclic_block;
    ci_max.p_cyclic_block = smz[ci_max_idx].p_cyclic_block;
    for (unsigned int c = 0; c < num_channels; c++) {
        sort(smz.begin(),
             smz.end(),
             [c](SequencingModel a, SequencingModel b) -> bool {
                 return a.channel_models[c]->p_bleach
                        < b.channel_models[c]->p_bleach;
             });
        ci_min.channel_models[c]->p_bleach =
                smz[ci_min_idx].channel_models[c]->p_bleach;
        ci_max.channel_models[c]->p_bleach =
                smz[ci_max_idx].channel_models[c]->p_bleach;
        sort(smz.begin(),
             smz.end(),
             [c](SequencingModel a, SequencingModel b) -> bool {
                 return a.channel_models[c]->p_dud < b.channel_models[c]->p_dud;
             });
        ci_min.channel_models[c]->p_dud =
                smz[ci_min_idx].channel_models[c]->p_dud;
        ci_max.channel_models[c]->p_dud =
                smz[ci_max_idx].channel_models[c]->p_dud;
    }
    cout << "lower bounds: " << ci_min.debug_string() << "\n";
    cout << "upper bounds: " << ci_max.debug_string() << "\n";
    cout << "worst bootstrap step size: " << worst_bootstrap_step_size << "\n";
    return best_log_l;
}

}  // namespace whatprot

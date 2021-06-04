/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "generate-dye-track.h"

// Standard C++ library headers:
#include <random>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

namespace {
using std::bernoulli_distribution;
using std::default_random_engine;
using std::vector;
}  // namespace

void generate_dye_track(const SequencingModel& seq_model,
                        const DyeSeq& dye_seq,
                        int num_timesteps,
                        int num_channels,
                        default_random_engine* generator,
                        DyeTrack* dye_track) {
    bernoulli_distribution edman_failure(seq_model.p_edman_failure);
    bernoulli_distribution detach_event(seq_model.p_detach);
    vector<bernoulli_distribution> dud_events;
    vector<bernoulli_distribution> bleach_events;
    for (int c = 0; c < num_channels; c++) {
        dud_events.emplace_back(seq_model.channel_models[c]->p_dud);
        bleach_events.emplace_back(seq_model.channel_models[c]->p_bleach);
    }
    DyeSeq ds(dye_seq);
    // Duds.
    for (int i = 0; i < ds.length; i++) {
        if (ds[i] != -1) {
            if (dud_events[ds[i]](*generator)) {
                ds[i] = -1;
            }
        }
    }
    // Track total counts for each channel.
    short* counts = new short[num_channels]();
    for (int i = 0; i < ds.length; i++) {
        if (ds[i] != -1) {
            counts[ds[i]]++;
        }
    }
    int e = 0;  // Successful Edman cycles.
    int t = 0;  // Timesteps.
    while (true) {
        // Record results.
        for (int c = 0; c < num_channels; c++) {
            (*dye_track)(t, c) = counts[c];
        }
        // Detach.
        if (detach_event(*generator)) {
            break;
        }
        t++;
        if (t >= num_timesteps) {
            break;
        }
        // Edman failures.
        if (!edman_failure(*generator)) {
            if (ds[e] != -1) {
                counts[ds[e]]--;
            }
            e++;
        }
        if (e >= ds.length) {
            break;
        }
        // Bleaching.
        for (int i = e; i < ds.length; i++) {
            if (ds[i] != -1) {
                if (bleach_events[ds[i]](*generator)) {
                    counts[ds[i]]--;
                    ds[i] = -1;
                }
            }
        }
    }
    while (t < num_timesteps) {
        for (int c = 0; c < num_channels; c++) {
            (*dye_track)(t, c) = 0;
        }
        t++;
    }
    delete[] counts;
}

}  // namespace whatprot

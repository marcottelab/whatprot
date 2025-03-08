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
                        unsigned int num_timesteps,
                        unsigned int num_channels,
                        default_random_engine* generator,
                        DyeTrack* dye_track) {
    bernoulli_distribution edman_failure(seq_model.p_edman_failure);
    vector<bernoulli_distribution> detach_events;
    for (unsigned int t = 0; t < num_timesteps; t++) {
        detach_events.emplace_back(seq_model.p_detach[t]);
    }
    bernoulli_distribution initial_block(seq_model.p_initial_block);
    bernoulli_distribution cyclic_block(seq_model.p_cyclic_block);
    vector<bernoulli_distribution> dud_events;
    vector<bernoulli_distribution> bleach_events;
    for (unsigned int c = 0; c < num_channels; c++) {
        dud_events.emplace_back(seq_model.channel_models[c]->p_dud);
        bleach_events.emplace_back(seq_model.channel_models[c]->p_bleach);
    }
    DyeSeq ds(dye_seq);
    // Duds.
    for (unsigned int i = 0; i < ds.length; i++) {
        if (ds[i] != -1) {
            if (dud_events[ds[i]](*generator)) {
                ds[i] = -1;
            }
        }
    }
    // Track total counts for each channel.
    short* counts = new short[num_channels]();
    for (unsigned int i = 0; i < ds.length; i++) {
        if (ds[i] != -1) {
            counts[ds[i]]++;
        }
    }
    unsigned int e = 0;  // Successful Edman cycles.
    bool block = false;
    if (initial_block(*generator)) {
        block = true;
    }
    for (unsigned int t = 0; t < num_timesteps; t++) {
        if (e < ds.length) {
            // Record results.
            for (unsigned int c = 0; c < num_channels; c++) {
                (*dye_track)(t, c) = counts[c];
            }
            // block n
            if (cyclic_block(*generator)) {
                block = true;
            }
            // Detach.
            if (detach_events[t](*generator)) {
                e = ds.length;
                continue;
            }
            // Edman failures.
            if (!block) {
                if (!edman_failure(*generator)) {
                    if (ds[e] != -1) {
                        counts[ds[e]]--;
                    }
                    e++;
                }
            }
            // Bleaching.
            for (unsigned int i = e; i < ds.length; i++) {
                if (ds[i] != -1) {
                    if (bleach_events[ds[i]](*generator)) {
                        counts[ds[i]]--;
                        ds[i] = -1;
                    }
                }
            }
        } else {
            // Record results.
            for (unsigned int c = 0; c < num_channels; c++) {
                (*dye_track)(t, c) = 0;
            }
        }
    }
    delete[] counts;
}

}  // namespace whatprot

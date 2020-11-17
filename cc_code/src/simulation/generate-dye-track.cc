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

// Local project headers:
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "common/error-model.h"

namespace fluoroseq {

namespace {
using std::bernoulli_distribution;
using std::default_random_engine;
}  // namespace

void generate_dye_track(const ErrorModel& error_model,
                        const DyeSeq& dye_seq,
                        int num_timesteps,
                        int num_channels,
                        default_random_engine* generator,
                        DyeTrack* dye_track) {
    bernoulli_distribution edman_failure(error_model.p_edman_failure);
    bernoulli_distribution dud_event(error_model.p_dud);
    bernoulli_distribution bleach_event(error_model.p_bleach);
    bernoulli_distribution detach_event(error_model.p_detach);
    DyeSeq ds(dye_seq);
    // Duds.
    for (int i = 0; i < ds.length; i++) {
        if (ds[i] != -1) {
            if (dud_event(*generator)) {
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
                if (bleach_event(*generator)) {
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

}  // namespace fluoroseq

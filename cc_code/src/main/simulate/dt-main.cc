/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "dt-main.h"

// Standard C++ library headers:
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/sourced-data.h"
#include "io/dye-seqs-io.h"
#include "io/dye-tracks-io.h"
#include "main/cmd-line-out.h"
#include "parameterization/model/sequencing-model.h"
#include "simulation/dedup-dye-tracks.h"
#include "simulation/generate-dye-tracks.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::atoi;
using std::default_random_engine;
using std::vector;
}  // namespace

int dt_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 7) {
        print_wrong_number_of_inputs();
        return EXIT_FAILURE;
    }
    unsigned int num_timesteps = atoi(argv[3]);
    int dye_tracks_per_peptide = atoi(argv[4]);
    char* dye_seqs_filename = argv[5];
    char* dye_tracks_filename = argv[6];

    double start_time;
    double end_time;

    start_time = wall_time();
    unsigned int num_channels;
    unsigned int total_num_dye_seqs;
    vector<SourcedData<DyeSeq, SourceCount<int>>> dye_seqs;
    read_dye_seqs(
            dye_seqs_filename, &num_channels, &total_num_dye_seqs, &dye_seqs);
    end_time = wall_time();
    print_read_dye_seqs(total_num_dye_seqs, end_time - start_time);

    start_time = wall_time();
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.06;
    seq_model.p_detach = 0.05;
    for (unsigned int c = 0; c < num_channels; c++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[c]->p_bleach = 0.05;
        seq_model.channel_models[c]->p_dud = 0.07;
        seq_model.channel_models[c]->bg_sigma = 0.00667;
        seq_model.channel_models[c]->mu = 1.0;
        seq_model.channel_models[c]->sigma = 0.16;
        seq_model.channel_models[c]->stuck_dye_ratio = 0.5;
        seq_model.channel_models[c]->p_stuck_dye_loss = 0.08;
    }
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    default_random_engine generator(time_based_seed());
    vector<SourcedData<DyeTrack, SourceCount<int>>> dye_tracks;
    generate_dye_tracks(seq_model,
                        dye_seqs,
                        num_timesteps,
                        num_channels,
                        dye_tracks_per_peptide,
                        &generator,
                        &dye_tracks);
    end_time = wall_time();
    print_finished_generating_dye_tracks(dye_tracks.size(),
                                         end_time - start_time);

    start_time = wall_time();
    vector<SourcedData<DyeTrack, SourceCountHitsList<int>>> deduped_dye_tracks;
    dedup_dye_tracks(
            num_timesteps, num_channels, &dye_tracks, &deduped_dye_tracks);
    end_time = wall_time();
    print_finished_deduping_dye_tracks(deduped_dye_tracks.size(),
                                       end_time - start_time);

    start_time = wall_time();
    write_dye_tracks(dye_tracks_filename,
                     num_timesteps,
                     num_channels,
                     deduped_dye_tracks);
    end_time = wall_time();
    print_finished_saving_results(end_time - start_time);

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);

    return 0;
}

}  // namespace whatprot

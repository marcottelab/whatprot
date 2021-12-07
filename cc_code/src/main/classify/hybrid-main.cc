/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "hybrid-main.h"

// Standard C++ library headers:
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// Local project headers:
#include "classifiers/hybrid-classifier.h"
#include "common/dye-track.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "io/dye-seqs-io.h"
#include "io/dye-tracks-io.h"
#include "io/radiometries-io.h"
#include "io/scored-classifications-io.h"
#include "main/cmd-line-out.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "util/delete.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::atof;
using std::atoi;
using std::string;
using std::vector;
}  // namespace

int hybrid_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 10) {
        print_wrong_number_of_inputs();
        return EXIT_FAILURE;
    }
    int k = atoi(argv[3]);
    double sig = atof(argv[4]);
    int h = atoi(argv[5]);
    char* dye_seqs_filename = argv[6];
    char* dye_tracks_filename = argv[7];
    char* radiometries_filename = argv[8];
    char* predictions_filename = argv[9];

    double start_time;
    double end_time;

    start_time = wall_time();
    unsigned int num_channels;
    unsigned int total_num_dye_seqs;  // redundant, not needed.
    vector<SourcedData<DyeSeq, SourceCount<int>>> dye_seqs;
    read_dye_seqs(
            dye_seqs_filename, &num_channels, &total_num_dye_seqs, &dye_seqs);
    end_time = wall_time();
    print_read_dye_seqs(dye_seqs.size(), end_time - start_time);

    start_time = wall_time();
    unsigned int num_timesteps;
    unsigned int duplicate_num_channels;  // also get this from dye seqs file
    vector<SourcedData<DyeTrack, SourceCountHitsList<int>>> dye_tracks;
    read_dye_tracks(dye_tracks_filename,
                    &num_timesteps,
                    &duplicate_num_channels,
                    &dye_tracks);
    end_time = wall_time();
    print_read_dye_tracks(dye_tracks.size(), end_time - start_time);

    start_time = wall_time();
    unsigned int duplicate_num_timesteps;  // also get this from dye track file.
    unsigned int triplicate_num_channels;  // see dye tracks and dye seqs files.
    unsigned int total_num_radiometries;  // num radiometries across all procs.
    vector<Radiometry> radiometries;
    read_radiometries(radiometries_filename,
                      &duplicate_num_timesteps,
                      &triplicate_num_channels,
                      &total_num_radiometries,
                      &radiometries);
    end_time = wall_time();
    print_read_radiometries(total_num_radiometries, end_time - start_time);

    start_time = wall_time();
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.06;
    seq_model.p_detach = 0.05;
    for (unsigned int c = 0; c < num_channels; c++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[c]->p_bleach = 0.05;
        seq_model.channel_models[c]->p_dud = 0.07;
        seq_model.channel_models[c]->bg_sig = 0.00667;
        seq_model.channel_models[c]->mu = 1.0;
        seq_model.channel_models[c]->sig = 0.16;
        seq_model.channel_models[c]->stuck_dye_ratio = 0.5;
        seq_model.channel_models[c]->p_stuck_dye_loss = 0.08;
    }
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = 3.0;
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    HybridClassifier classifier(num_timesteps,
                                num_channels,
                                seq_model,
                                seq_settings,
                                k,
                                sig,
                                &dye_tracks,
                                h,
                                dye_seqs);
    end_time = wall_time();
    print_built_classifier(end_time - start_time);

    start_time = wall_time();
    vector<ScoredClassification> results = classifier.classify(radiometries);
    end_time = wall_time();
    print_finished_classification(end_time - start_time);

    start_time = wall_time();
    write_scored_classifications(
            predictions_filename, total_num_radiometries, results);
    end_time = wall_time();
    print_finished_saving_results(end_time - start_time);

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);

    return 0;
}

}  // namespace whatprot

// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#include "hybrid_main.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "classifiers/hybrid_classifier.h"
#include "common/approximation_model.h"
#include "common/dye_track.h"
#include "common/error_model.h"
#include "common/radiometry.h"
#include "common/scored_classification.h"
#include "common/sourced_data.h"
#include "io/dye_seqs_io.h"
#include "io/dye_tracks_io.h"
#include "io/radiometries_io.h"
#include "io/scored_classifications_io.h"
#include "main/cmd_line_out.h"
#include "util/delete.h"
#include "util/time.h"

namespace fluoroseq {

namespace {
using std::atoi;
using std::string;
using std::vector;
}  // namespace

int hybrid_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 7) {
        print_wrong_number_of_inputs();
        return EXIT_FAILURE;
    }
    int h = atoi(argv[2]);
    char* dye_seqs_filename = argv[3];
    char* dye_tracks_filename = argv[4];
    char* radiometries_filename = argv[5];
    char* predictions_filename = argv[6];

    double start_time;
    double end_time;

    start_time = wall_time();
    ErrorModel error_model(.06,  // p_edman_failure
                           .05,  // p_detach
                           .05,  // p_bleach
                           .07,  // p_dud
                           DistributionType::LOGNORMAL,
                           1.0,  // mu
                           .16);  // sigma
    ApproximationModel approximation_model(16);
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    int num_channels;
    int total_num_dye_seqs;  // redundant, not needed.
    vector<SourcedData<DyeSeq, SourceCount<int>>> dye_seqs;
    read_dye_seqs(dye_seqs_filename,
                  &num_channels,
                  &total_num_dye_seqs,
                  &dye_seqs,
                  MpiReadMode::BROADCAST);
    end_time = wall_time();
    print_read_dye_seqs(dye_seqs.size(), end_time - start_time);

    start_time = wall_time();
    int num_timesteps;
    int duplicate_num_channels;  // also get this from dye seqs file
    int num_dye_tracks;
    vector<SourcedData<DyeTrack, SourceCountHitsList<int>>> dye_tracks;
    read_dye_tracks(dye_tracks_filename,
                    &num_timesteps,
                    &duplicate_num_channels,
                    &dye_tracks);
    end_time = wall_time();
    print_read_dye_tracks(num_dye_tracks, end_time - start_time);

    start_time = wall_time();
    int duplicate_num_timesteps;  // also get this from dye track file.
    int triplicate_num_channels;  // also in dye tracks and dye seqs files.
    int total_num_radiometries;  // number of radiometries across all procs.
    vector<Radiometry> radiometries;
    read_radiometries(radiometries_filename,
                      &duplicate_num_timesteps,
                      &triplicate_num_channels,
                      &total_num_radiometries,
                      &radiometries);
    end_time = wall_time();
    print_read_radiometries(total_num_radiometries, end_time - start_time);

    start_time = wall_time();
    HybridClassifier classifier(num_timesteps,
                                num_channels,
                                error_model,
                                approximation_model,
                                10,  // k
                                dye_tracks,
                                h,
                                dye_seqs);
    end_time = wall_time();
    print_built_classifier(end_time - start_time);

    start_time = wall_time();
    vector<ScoredClassification> results = classifier.classify(radiometries);
    end_time = wall_time();
    print_finished_classification(end_time - start_time);
    
    start_time = wall_time();
    write_scored_classifications(predictions_filename,
                                 total_num_radiometries,
                                 results);
    end_time = wall_time();
    print_finished_saving_results(end_time - start_time);

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);

    return 0;
}

}  // namespace fluoroseq

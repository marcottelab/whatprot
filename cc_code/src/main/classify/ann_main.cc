// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#include "ann_main.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "classifiers/kwann_classifier.h"
#include "common/dye_track.h"
#include "common/error_model.h"
#include "common/radiometry.h"
#include "common/scored_classification.h"
#include "common/sourced_data.h"
#include "io/dye_tracks_io.h"
#include "io/radiometries_io.h"
#include "io/scored_classifications_io.h"
#include "main/cmd_line_out.h"
#include "util/time.h"

namespace fluoroseq {

namespace {
using std::atoi;
using std::vector;
}  // namespace

int ann_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 6) {
        print_wrong_number_of_inputs();
        return EXIT_FAILURE;
    }
    int k = atoi(argv[2]);
    char* dye_tracks_filename = argv[3];
    char* radiometries_filename = argv[4];
    char* predictions_filename = argv[5];

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
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    int num_timesteps;
    int num_channels;
    vector<SourcedData<DyeTrack, SourceCountHitsList<int>>> dye_tracks;
    read_dye_tracks(dye_tracks_filename,
                    &num_timesteps,
                    &num_channels,
                    &dye_tracks);
    end_time = wall_time();
    print_read_dye_tracks(dye_tracks.size(), end_time - start_time);

    start_time = wall_time();
    int duplicate_num_timesteps;  // also get this from dye track file.
    int duplicate_num_channels;  // also get this from dye track file.
    int total_num_radiometries;  // number of radiometries across all procs.
    vector<Radiometry> radiometries;
    read_radiometries(radiometries_filename,
                      &duplicate_num_timesteps,
                      &duplicate_num_channels,
                      &total_num_radiometries,
                      &radiometries);
    end_time = wall_time();
    print_read_radiometries(total_num_radiometries, end_time - start_time);

    start_time = wall_time();
    KWANNClassifier classifier(num_timesteps,
                               num_channels,
                               error_model,
                               k,
                               dye_tracks);
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

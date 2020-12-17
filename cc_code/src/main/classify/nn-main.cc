/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// For MPI version, define compiler macro USE_MPI when building.

// Defining symbols from header:
#include "nn-main.h"

// Standard C++ library headers:
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// Local project headers:
#include "classifiers/nn-classifier.h"
#include "common/dye-track.h"
#include "common/error-model.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "io/dye-tracks-io.h"
#include "io/radiometries-io.h"
#include "io/scored-classifications-io.h"
#include "main/cmd-line-out.h"
#include "util/time.h"

namespace fluoroseq {

namespace {
using std::atof;
using std::atoi;
using std::vector;
}  // namespace

int nn_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 8) {
        print_wrong_number_of_inputs();
        return EXIT_FAILURE;
    }
    int k = atoi(argv[3]);
    double sigma = atof(argv[4]);
    char* dye_tracks_filename = argv[5];
    char* radiometries_filename = argv[6];
    char* predictions_filename = argv[7];

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
    read_dye_tracks(
            dye_tracks_filename, &num_timesteps, &num_channels, &dye_tracks);
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
    NNClassifier classifier(
            num_timesteps, num_channels, error_model, k, sigma, &dye_tracks);
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

}  // namespace fluoroseq

/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "hmm-main.h"

// Standard C++ library headers:
#include <iostream>
#include <string>
#include <vector>

// Local project headers:
#include "classifiers/hmm-classifier.h"
#include "common/dye-seq.h"
#include "common/error-model.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "io/dye-seqs-io.h"
#include "io/radiometries-io.h"
#include "io/scored-classifications-io.h"
#include "main/cmd-line-out.h"
#include "util/delete.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::string;
using std::vector;
}  // namespace

int hmm_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 6) {
        print_wrong_number_of_inputs();
        return EXIT_FAILURE;
    }
    char* dye_seqs_filename = argv[3];
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
                           0.0,  // mu
                           .16,  // sigma
                           0.5);  // stuck_dye_ratio
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    int num_channels;
    int total_num_dye_seqs;  // redundant, not needed.
    vector<SourcedData<DyeSeq, SourceCount<int>>> dye_seqs;
    read_dye_seqs(
            dye_seqs_filename, &num_channels, &total_num_dye_seqs, &dye_seqs);
    end_time = wall_time();
    print_read_dye_seqs(dye_seqs.size(), end_time - start_time);

    start_time = wall_time();
    int num_timesteps;
    int duplicate_num_channels;  // also get this from dye seq file.
    int total_num_radiometries;  // number of radiometries across all procs.
    vector<Radiometry> radiometries;
    read_radiometries(radiometries_filename,
                      &num_timesteps,
                      &duplicate_num_channels,
                      &total_num_radiometries,
                      &radiometries);
    end_time = wall_time();
    print_read_radiometries(total_num_radiometries, end_time - start_time);

    start_time = wall_time();
    HMMClassifier classifier(
            num_timesteps, num_channels, error_model, dye_seqs);
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

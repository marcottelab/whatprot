// Author: Matthew Beauregard Smith (UT Austin)
// Simple application to read in a TSV file of dye seqs, a TSV file of
// radiometries, and write predicted classifications to a TSV file.
//
// For MPI version, define compiler macro USE_MPI when building.
#include "hmm_main.h"

#include <iostream>
#include <string>

#include "classifiers/fwd_alg_classifier.h"
#include "common/approximation_model.h"
#include "common/dye_seq.h"
#include "common/error_model.h"
#include "common/radiometry.h"
#include "common/scored_classification.h"
#include "common/sourced_data.h"
#include "io/dye_seqs_io.h"
#include "io/radiometries_io.h"
#include "io/scored_classifications_io.h"
#include "main/cmd_line_out.h"
#include "util/time.h"

namespace fluoroseq {

namespace {
using std::string;
}

int hmm_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 5) {
        print_wrong_number_of_inputs();
        return EXIT_FAILURE;
    }
    char* dye_seqs_filename = argv[2];
    char* radiometries_filename = argv[3];
    char* predictions_filename = argv[4];

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
    int num_dye_seqs;
    SourcedData<DyeSeq*, SourceCount<int>*>** dye_seqs;
    read_dye_seqs(dye_seqs_filename,
                  &num_channels,
                  &num_dye_seqs,
                  &dye_seqs);
    end_time = wall_time();
    print_read_dye_seqs(num_dye_seqs, end_time - start_time);

    start_time = wall_time();
    int num_timesteps;
    int duplicate_num_channels;  // also get this from dye seq file.
    int total_num_radiometries;  // number of radiometries across all procs.
    int num_radiometries;  // number of radiometries in this proc.
    Radiometry** radiometries;
    read_radiometries(radiometries_filename,
                      &num_timesteps,
                      &duplicate_num_channels,
                      &total_num_radiometries,
                      &num_radiometries,
                      &radiometries);
    end_time = wall_time();
    print_read_radiometries(total_num_radiometries, end_time - start_time);

    start_time = wall_time();
    FwdAlgClassifier classifier(num_timesteps,
                                num_channels,
                                error_model,
                                approximation_model,
                                num_dye_seqs,
                                dye_seqs);
    end_time = wall_time();
    print_built_classifier(end_time - start_time);

    start_time = wall_time();
    ScoredClassification* results = classifier.classify(num_radiometries,
                                                         radiometries);
    end_time = wall_time();
    print_finished_classification(end_time - start_time);
    
    start_time = wall_time();
    write_scored_classifications(predictions_filename,
                                 total_num_radiometries,
                                 num_radiometries,
                                 results);
    end_time = wall_time();
    print_finished_saving_results(end_time - start_time);

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);

    return 0;
}

}  // namespace fluoroseq

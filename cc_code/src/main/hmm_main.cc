// Author: Matthew Beauregard Smith (UT Austin)
// Simple application to read in a TSV file of dye seqs, a TSV file of
// radiometries, and write predicted classifications to a TSV file.
#include "hmm_main.h"

#include <iostream>
#include <string>

#include "common/approximation_model.h"
#include "common/dye_seq.h"
#include "common/error_model.h"
#include "common/radiometry.h"
#include "common/scored_classification.h"
#include "common/sourced_data.h"
#include "io/dye_seqs_io.h"
#include "io/radiometries_io.h"
#include "io/scored_classifications_io.h"
#include "util/time.h"

#ifdef _OPENMP
#include "classifiers/omp_fwd_alg_classifier.h"
#else
#include "classifiers/fwd_alg_classifier.h"
#endif

namespace fluoroseq {

namespace {
using std::cout;
using std::string;
}

int hmm_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 5) {
        cout << "wrong number of inputs\n";
        return EXIT_FAILURE;
    }
    char* dye_seqs_filename = argv[2];
    char* radiometries_filename = argv[3];
    char* predictions_filename = argv[4];

    double start_time;
    double end_time;

    #ifdef _OPENMP
    cout << "Using OpenMP\n";
    cout << "    Number of threads: " << omp_get_max_threads() << "\n";
    #endif

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
    cout << "Built error model and approximation model.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";

    start_time = wall_time();
    int num_channels;
    int num_dye_seqs;
    SourcedData<DyeSeq*, SourceWithCount<int>*>** dye_seqs;
    read_dye_seqs(dye_seqs_filename,
                  &num_channels,
                  &num_dye_seqs,
                  &dye_seqs);
    end_time = wall_time();
    cout << "Read from dye seqs file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Number of dye seqs: " << num_dye_seqs << "\n";

    start_time = wall_time();
    int num_timesteps;
    int duplicate_num_channels;  // also get this from dye seq file.
    int num_radiometries;
    Radiometry** radiometries;
    read_radiometries(radiometries_filename,
                      &num_timesteps,
                      &duplicate_num_channels,
                      &num_radiometries,
                      &radiometries);
    end_time = wall_time();
    cout << "Read from radiometries file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Number of radiometries: " << num_radiometries << "\n";

    start_time = wall_time();
    #ifdef _OPENMP
    OMPFwdAlgClassifier classifier(num_timesteps,
                                   num_channels,
                                   error_model,
                                   approximation_model,
                                   num_dye_seqs,
                                   dye_seqs);
    #else
    FwdAlgClassifier classifier(num_timesteps,
                                num_channels,
                                error_model,
                                approximation_model,
                                num_dye_seqs,
                                dye_seqs);
    #endif
    end_time = wall_time();
    cout << "Constructed classifier.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";

    start_time = wall_time();
    ScoredClassification* results = classifier.classify(num_radiometries,
                                                         radiometries);
    end_time = wall_time();
    cout << "Got results.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Total number of runs: "
         << num_dye_seqs * num_radiometries
         << "\n";
    cout << "    Time per run in seconds: "
         << (end_time - start_time) / (num_dye_seqs * num_radiometries)
         << "\n";
    
    start_time = wall_time();
    write_scored_classifications(predictions_filename,
                                 num_radiometries,
                                 results);
    end_time = wall_time();
    cout << "Wrote to file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";

    double total_end_time = wall_time();
    cout << "Totals:\n";
    cout << "    Time in seconds: "
         << total_end_time - total_start_time
         << "\n";
    cout << "    Total number of runs: "
         << num_dye_seqs * num_radiometries
         << "\n";
    cout << "    Time per run in seconds: "
         << (total_end_time - total_start_time)
            / (num_dye_seqs * num_radiometries)
         << "\n";

    return 0;
}

}  // namespace fluoroseq

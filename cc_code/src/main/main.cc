// Author: Matthew Beauregard Smith (UT Austin)
// Simple application to read in a TSV file of dye seqs, a TSV file of
// radiometries, and write predicted classifications to a TSV file.

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

#ifdef _OPENMP
#include <omp.h>
#include "classifiers/omp_fwd_alg_classifier.h"
#else
#include <ctime>
#include "classifiers/fwd_alg_classifier.h"
#endif

namespace {
using fluoroseq::read_dye_seqs;
using fluoroseq::read_radiometries;
using fluoroseq::write_scored_classifications;
using fluoroseq::ApproximationModel;
using fluoroseq::DistributionType;
using fluoroseq::DyeSeq;
using fluoroseq::ErrorModel;
using fluoroseq::Radiometry;
using fluoroseq::ScoredClassification;
using fluoroseq::SourcedData;
using fluoroseq::SourceWithCount;
using std::cout;
using std::string;
#ifdef _OPENMP
using fluoroseq::OMPFwdAlgClassifier;
#else
using fluoroseq::FwdAlgClassifier;
#endif
}

double wtime();

int main(int argc, char** argv) {
    double total_start_time = wtime();

    if (argc != 4) {
        cout << "wrong number of inputs\n";
        return EXIT_FAILURE;
    }
    char* dye_seqs_filename = argv[1];
    char* radiometries_filename = argv[2];
    char* predictions_filename = argv[3];

    double start_time;
    double end_time;

    #ifdef _OPENMP
    cout << "Using OpenMP\n";
    cout << "    Number of threads: " << omp_get_max_threads() << "\n";
    #endif

    start_time = wtime();
    ErrorModel error_model(.06,  // p_edman_failure
                           .05,  // p_detach
                           .05,  // p_bleach
                           .07,  // p_dud
                           DistributionType::LOGNORMAL,
                           1.0,  // mu
                           .16);  // sigma
    ApproximationModel approximation_model(0);
    end_time = wtime();
    cout << "Built error model and approximation model.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";

    start_time = wtime();
    int num_channels;
    int num_dye_seqs;
    SourcedData<DyeSeq*, SourceWithCount<int>>** dye_seqs;
    read_dye_seqs(dye_seqs_filename,
                  &num_channels,
                  &num_dye_seqs,
                  &dye_seqs);
    end_time = wtime();
    cout << "Read from dye seqs file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Number of dye seqs: " << num_dye_seqs << "\n";

    start_time = wtime();
    int num_timesteps;
    int duplicate_num_channels;  // also get this from dye seq file.
    int num_radiometries;
    Radiometry** radiometries = new Radiometry*[num_radiometries];
    read_radiometries(radiometries_filename,
                      &num_timesteps,
                      &duplicate_num_channels,
                      &num_radiometries,
                      &radiometries);
    end_time = wtime();
    cout << "Read from radiometries file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Number of radiometries: " << num_radiometries << "\n";

    start_time = wtime();
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
    end_time = wtime();
    cout << "Constructed classifier.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";

    start_time = wtime();
    ScoredClassification* results = classifier.classify(num_radiometries,
                                                         radiometries);
    end_time = wtime();
    cout << "Got results.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Total number of runs: "
         << num_dye_seqs * num_radiometries
         << "\n";
    cout << "    Time per run in seconds: "
         << (end_time - start_time) / (num_dye_seqs * num_radiometries)
         << "\n";
    
    start_time = wtime();
    write_scored_classifications(predictions_filename,
                                 num_radiometries,
                                 results);
    end_time = wtime();
    cout << "Wrote to file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";

    double total_end_time = wtime();
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

double wtime() {
    #ifdef _OPENMP
    return omp_get_wtime();
    #else
    return (double) clock() / (double) CLOCKS_PER_SEC;
    #endif
}

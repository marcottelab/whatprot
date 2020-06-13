// Author: Matthew Beauregard Smith (UT Austin)
#include "hybrid_main.h"

#include <iostream>
#include <string>

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
#include "util/time.h"

namespace fluoroseq {

namespace {
using std::cout;
using std::string;
}  // namespace

int hybrid_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 6) {
        cout << "wrong number of inputs\n";
        return EXIT_FAILURE;
    }
    char* dye_seqs_filename = argv[2];
    char* dye_tracks_filename = argv[3];
    char* radiometries_filename = argv[4];
    char* predictions_filename = argv[5];

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
    int duplicate_num_channels;  // also get this from dye seqs file
    int num_dye_tracks;
    SourcedData<DyeTrack*, SourceCountMap<int>*>** dye_tracks;
    read_dye_tracks(dye_tracks_filename,
                    &num_timesteps,
                    &duplicate_num_channels,
                    &num_dye_tracks,
                    &dye_tracks);
    end_time = wall_time();
    cout << "Read from dye tracks file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Number of dye seqs: " << num_dye_tracks << "\n";

    start_time = wall_time();
    int duplicate_num_timesteps;  // also get this from dye track file.
    int triplicate_num_channels;  // also in dye tracks and dye seqs files.
    int num_radiometries;
    Radiometry** radiometries;
    read_radiometries(radiometries_filename,
                      &duplicate_num_timesteps,
                      &triplicate_num_channels,
                      &num_radiometries,
                      &radiometries);
    end_time = wall_time();
    cout << "Read from radiometries file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Number of radiometries: " << num_radiometries << "\n";

    start_time = wall_time();
    HybridClassifier classifier(num_timesteps,
                                num_channels,
                                error_model,
                                approximation_model,
                                1000,  // k
                                num_dye_tracks,
                                dye_tracks,
                                100,  // h
                                num_dye_seqs,
                                dye_seqs);
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
         << num_dye_tracks * num_radiometries
         << "\n";
    cout << "    Time per run in seconds: "
         << (end_time - start_time) / (num_dye_tracks * num_radiometries)
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
         << num_dye_tracks * num_radiometries
         << "\n";
    cout << "    Time per run in seconds: "
         << (total_end_time - total_start_time)
            / (num_dye_tracks * num_radiometries)
         << "\n";

    return 0;
}

}  // namespace fluoroseq

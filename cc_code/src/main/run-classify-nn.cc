/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "run-classify-nn.h"

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
#include "classifiers/nn-classifier.h"
#include "common/dye-track.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "io/dye-tracks-io.h"
#include "io/radiometries-io.h"
#include "io/scored-classifications-io.h"
#include "main/cmd-line-out.h"
#include "parameterization/model/sequencing-model.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::string;
using std::vector;
}  // namespace

void run_classify_nn(int k,
                     double sig,
                     string dye_tracks_filename,
                     string radiometries_filename,
                     string predictions_filename) {
    double total_start_time = wall_time();

    double start_time;
    double end_time;

    start_time = wall_time();
    unsigned int num_timesteps;
    unsigned int num_channels;
    vector<SourcedData<DyeTrack, SourceCountHitsList<int>>> dye_tracks;
    read_dye_tracks(
            dye_tracks_filename, &num_timesteps, &num_channels, &dye_tracks);
    end_time = wall_time();
    print_read_dye_tracks(dye_tracks.size(), end_time - start_time);

    start_time = wall_time();
    unsigned int duplicate_num_timesteps;  // also get this from dye track file.
    unsigned int duplicate_num_channels;  // also get this from dye track file.
    unsigned int total_num_radiometries;  // num radiometries across all procs.
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
            num_timesteps, num_channels, k, sig, &dye_tracks);
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
}

}  // namespace whatprot

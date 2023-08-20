/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "run-score-hmm.h"

// Standard C++ library headers:
#include <limits>
#include <string>
#include <vector>

// Local project headers:
#include "classifiers/hmm-classifier.h"
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "common/sourced-data.h"
#include "io/dye-seqs-io.h"
#include "io/radiometries-io.h"
#include "io/scores-io.h"
#include "main/cmd-line-out.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "util/delete.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::string;
using std::vector;
}  // namespace

void run_score_hmm(double hmm_pruning_cutoff,
                   string dye_seqs_filename,
                   string radiometries_filename,
                   string scores_filename) {
    double total_start_time = wall_time();

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
    unsigned int duplicate_num_channels;  // also get this from dye seq file.
    unsigned int total_num_radiometries;  // num radiometries across all procs.
    vector<Radiometry> radiometries;
    read_radiometries(radiometries_filename,
                      &num_timesteps,
                      &duplicate_num_channels,
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
    seq_settings.dist_cutoff = hmm_pruning_cutoff;
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    HMMClassifier classifier(
            num_timesteps, num_channels, seq_model, seq_settings, dye_seqs);
    end_time = wall_time();
    print_built_classifier(end_time - start_time);

    start_time = wall_time();
    vector<vector<double>> all_scores = classifier.score(radiometries);
    end_time = wall_time();
    print_finished_classification(end_time - start_time);

    start_time = wall_time();
    write_scores(scores_filename, all_scores);
    end_time = wall_time();
    print_finished_saving_results(end_time - start_time);

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);
}

}  // namespace whatprot

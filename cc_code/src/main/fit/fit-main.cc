/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Standard C++ library headers:
#include <cstdlib>
#include <cstring>

// Local project headers:
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "fitters/hmm-fitter.h"
#include "io/radiometries-io.h"
#include "main/cmd-line-out.h"
#include "parameterization/model/sequencing-model.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::atof;
using std::atoi;
using std::vector;
}  // namespace

int fit_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc < 5) {
        print_bad_inputs();
        return 1;
    }
    double stopping_threshold = atof(argv[2]);
    char* dye_seq_string = argv[3];
    char* radiometries_filename = argv[4];

    double start_time;
    double end_time;

    start_time = wall_time();
    int num_timesteps;
    int num_channels;
    int total_num_radiometries;
    vector<Radiometry> radiometries;
    read_radiometries(radiometries_filename,
                      &num_timesteps,
                      &num_channels,
                      &total_num_radiometries,
                      &radiometries);
    end_time = wall_time();
    print_read_radiometries(total_num_radiometries, end_time - start_time);

    start_time = wall_time();
    DyeSeq dye_seq(num_channels, dye_seq_string);
    end_time = wall_time();

    start_time = wall_time();
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.06;
    seq_model.p_detach = 0.05;
    for (int c = 0; c < num_channels; c++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[c]->p_bleach = 0.05;
        seq_model.channel_models[c]->p_dud = 0.07;
        seq_model.channel_models[c]->mu = 0.0;
        seq_model.channel_models[c]->sigma = 0.16;
        seq_model.channel_models[c]->stuck_dye_ratio = 0.5;
        seq_model.channel_models[c]->p_stuck_dye_loss = 0.08;
    }
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    HMMFitter fitter(num_timesteps,
                     num_channels,
                     stopping_threshold,
                     seq_model,
                     dye_seq);
    fitter.fit(radiometries);
    end_time = wall_time();

    return 0;
}

}  // namespace whatprot

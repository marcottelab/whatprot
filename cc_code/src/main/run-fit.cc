/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "run-fit.h"

// Standard C++ library headers:
#include <limits>
#include <string>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "fitters/hmm-fitter.h"
#include "io/radiometries-io.h"
#include "main/cmd-line-out.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::string;
using std::vector;
}  // namespace

void run_fit(double stopping_threshold,
             string dye_seq_string,
             std::string seq_params_filename,
             string radiometries_filename) {
    double total_start_time = wall_time();

    double start_time;
    double end_time;

    start_time = wall_time();
    unsigned int num_timesteps;
    unsigned int num_channels;
    unsigned int total_num_radiometries;
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
    SequencingModel seq_model(seq_params_filename);
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    HMMFitter fitter(num_timesteps,
                     num_channels,
                     stopping_threshold,
                     seq_model,
                     seq_settings,
                     dye_seq);
    fitter.fit(radiometries);
    print_finished_parameter_fitting(end_time - start_time);
    end_time = wall_time();

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);
}

}  // namespace whatprot

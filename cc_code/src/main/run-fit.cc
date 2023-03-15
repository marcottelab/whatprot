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
#include "fitters/bootstrap-fit.h"
#include "fitters/hmm-fitter.h"
#include "io/params-io.h"
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
             string seq_params_filename,
             string radiometries_filename,
             unsigned int num_bootstrap,
             double confidence_interval,
             string results_filename) {
    double total_start_time = wall_time();

    double start_time;
    double end_time;

    start_time = wall_time();
    SequencingModel true_seq_model(seq_params_filename);
    SequencingModel seq_model = true_seq_model.with_mu_as_one();
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    unsigned int num_timesteps;
    unsigned int num_channels;
    unsigned int total_num_radiometries;
    vector<Radiometry> radiometries;
    read_radiometries(radiometries_filename,
                      true_seq_model,
                      &num_timesteps,
                      &num_channels,
                      &total_num_radiometries,
                      &radiometries);
    end_time = wall_time();
    print_read_radiometries(total_num_radiometries, end_time - start_time);

    start_time = wall_time();
    DyeSeq dye_seq(num_channels, dye_seq_string);
    end_time = wall_time();

    if (confidence_interval == 0.0) {
        start_time = wall_time();
        HMMFitter fitter(num_timesteps,
                         num_channels,
                         stopping_threshold,
                         seq_model,
                         seq_settings,
                         dye_seq);
        SequencingModel fitted_seq_model;
        double log_l = fitter.fit(radiometries, &fitted_seq_model);
        end_time = wall_time();
        print_finished_parameter_fitting(end_time - start_time);

        print_parameter_results(fitted_seq_model, log_l);
    } else {
        start_time = wall_time();
        vector<SequencingModel> seq_models;
        vector<double> log_ls;
        bootstrap_fit(num_timesteps,
                      num_channels,
                      stopping_threshold,
                      seq_model,
                      seq_settings,
                      dye_seq,
                      radiometries,
                      num_bootstrap,
                      confidence_interval,
                      &seq_models,
                      &log_ls);
        end_time = wall_time();
        print_finished_parameter_fitting(end_time - start_time);

        if (results_filename != "") {
            start_time = wall_time();
            write_params(results_filename, num_channels, seq_models, log_ls);
            end_time = wall_time();
            print_finished_saving_results(end_time - start_time);
        }
    }

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);
}

}  // namespace whatprot

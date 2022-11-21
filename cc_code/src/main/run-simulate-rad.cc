/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "run-simulate-rad.h"

// Standard C++ library headers:
#include <random>
#include <string>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/sourced-data.h"
#include "io/dye-seqs-io.h"
#include "io/radiometries-io.h"
#include "main/cmd-line-out.h"
#include "parameterization/model/sequencing-model.h"
#include "simulation/generate-radiometries.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::default_random_engine;
using std::string;
using std::vector;
}  // namespace

void run_simulate_rad(unsigned int num_timesteps,
                      unsigned int num_to_generate,
                      string seq_params_filename,
                      string dye_seqs_filename,
                      string radiometries_filename,
                      string ys_filename) {
    double total_start_time = wall_time();

    double start_time;
    double end_time;

    start_time = wall_time();
    SequencingModel true_seq_model(seq_params_filename);
    SequencingModel seq_model = true_seq_model.with_mu_as_one();
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    unsigned int num_channels;
    unsigned int total_num_dye_seqs;
    vector<SourcedData<DyeSeq, SourceCount<int>>> dye_seqs;
    read_dye_seqs(
            dye_seqs_filename, &num_channels, &total_num_dye_seqs, &dye_seqs);
    end_time = wall_time();
    print_read_dye_seqs(total_num_dye_seqs, end_time - start_time);

    start_time = wall_time();
    default_random_engine generator(time_based_seed());
    vector<SourcedData<Radiometry, SourceCount<int>>> radiometries;
    generate_radiometries(seq_model,
                          dye_seqs,
                          num_timesteps,
                          num_channels,
                          num_to_generate,
                          &generator,
                          &radiometries);
    end_time = wall_time();
    print_finished_generating_radiometries(radiometries.size(),
                                           end_time - start_time);

    start_time = wall_time();
    write_radiometries(radiometries_filename,
                       true_seq_model,
                       num_timesteps,
                       num_channels,
                       radiometries);
    write_ys(ys_filename, radiometries);
    end_time = wall_time();
    print_finished_saving_results(end_time - start_time);

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);
}

}  // namespace whatprot

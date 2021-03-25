/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "rad-main.h"

// Standard C++ library headers:
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/error-model.h"
#include "common/sourced-data.h"
#include "io/dye-seqs-io.h"
#include "io/radiometries-io.h"
#include "main/cmd-line-out.h"
#include "simulation/generate-radiometries.h"
#include "util/time.h"

namespace whatprot {

namespace {
using std::atoi;
using std::default_random_engine;
using std::vector;
}  // namespace

int rad_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 8) {
        print_wrong_number_of_inputs();
        return EXIT_FAILURE;
    }
    int num_timesteps = atoi(argv[3]);
    int radiometries_per_peptide = atoi(argv[4]);
    char* dye_seqs_filename = argv[5];
    char* radiometries_filename = argv[6];
    char* ys_filename = argv[7];

    double start_time;
    double end_time;

    start_time = wall_time();
    ErrorModel error_model(.06,  // p_edman_failure
                           .05,  // p_detach
                           .05,  // p_bleach
                           .07,  // p_dud
                           DistributionType::LOGNORMAL,
                           0.0,  // mu
                           .16,  // sigma
                           0.5,  // stuck_dye_ratio
                           .08);  // p_stuck_dye_loss
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    int num_channels;
    int total_num_dye_seqs;
    vector<SourcedData<DyeSeq, SourceCount<int>>> dye_seqs;
    read_dye_seqs(
            dye_seqs_filename, &num_channels, &total_num_dye_seqs, &dye_seqs);
    end_time = wall_time();
    print_read_dye_seqs(total_num_dye_seqs, end_time - start_time);

    start_time = wall_time();
    default_random_engine generator(time_based_seed());
    vector<SourcedData<Radiometry, SourceCount<int>>> radiometries;
    generate_radiometries(error_model,
                          dye_seqs,
                          num_timesteps,
                          num_channels,
                          radiometries_per_peptide,
                          &generator,
                          &radiometries);
    end_time = wall_time();
    print_finished_generating_radiometries(radiometries.size(),
                                           end_time - start_time);

    start_time = wall_time();
    write_radiometries(
            radiometries_filename, num_timesteps, num_channels, radiometries);
    write_ys(ys_filename, radiometries);
    end_time = wall_time();
    print_finished_saving_results(end_time - start_time);

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);

    return 0;
}

}  // namespace whatprot

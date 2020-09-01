// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#include "rad_main.h"

#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "common/dye_seq.h"
#include "common/error_model.h"
#include "common/sourced_data.h"
#include "io/dye_seqs_io.h"
#include "io/radiometries_io.h"
#include "main/cmd_line_out.h"
#include "simulation/generate_radiometries.h"
#include "util/time.h"

namespace fluoroseq {

namespace {
using std::atoi;
using std::default_random_engine;
using std::vector;
}  // namespace

int rad_main(int argc, char** argv) {
    double total_start_time = wall_time();

    if (argc != 6) {
        print_wrong_number_of_inputs();
        return EXIT_FAILURE;
    }
    int num_timesteps = atoi(argv[2]);
    int radiometries_per_dye_seq = atoi(argv[3]);
    char* dye_seqs_filename = argv[4];
    char* radiometries_filename = argv[5];

    double start_time;
    double end_time;

    start_time = wall_time();
    ErrorModel error_model(.06,  // p_edman_failure
                           .05,  // p_detach
                           .05,  // p_bleach
                           .07,  // p_dud
                           DistributionType::LOGNORMAL,
                           1.0,  // mu
                           .16);  // sigma
    end_time = wall_time();
    print_finished_basic_setup(end_time - start_time);

    start_time = wall_time();
    int num_channels;
    int total_num_dye_seqs;
    vector<SourcedData<DyeSeq, SourceCount<int>>> dye_seqs;
    read_dye_seqs(dye_seqs_filename,
                  &num_channels,
                  &total_num_dye_seqs,
                  &dye_seqs,
                  MpiReadMode::SCATTER);
    end_time = wall_time();
    print_read_dye_seqs(total_num_dye_seqs, end_time - start_time);

    start_time = wall_time();
    default_random_engine generator(time_based_seed());
    vector<SourcedData<Radiometry, SourceCount<int>>> radiometries;
    generate_radiometries(error_model,
                          dye_seqs,
                          num_timesteps,
                          num_channels,
                          radiometries_per_dye_seq,
                          &generator,
                          &radiometries);
    end_time = wall_time();
    print_finished_generating_radiometries(
            total_num_dye_seqs * radiometries_per_dye_seq,
            end_time - start_time);

    start_time = wall_time();
    write_radiometries(radiometries_filename,
                       num_timesteps,
                       num_channels,
                       total_num_dye_seqs,  // num groups
                       radiometries_per_dye_seq,  // group size
                       radiometries);
    end_time = wall_time();
    print_finished_saving_results(end_time - start_time);

    double total_end_time = wall_time();
    print_total_time(total_end_time - total_start_time);

    return 0;
}

}  // namespace fluoroseq

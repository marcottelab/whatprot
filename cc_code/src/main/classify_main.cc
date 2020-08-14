// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.

#include <cstring>

#ifdef USE_MPI
#include <mpi.h>
#endif  // USE_MPI

#include "main/ann_main.h"
#include "main/cmd_line_out.h"
#include "main/hmm_main.h"
#include "main/hybrid_main.h"

namespace {
using fluoroseq::ann_main;
using fluoroseq::hmm_main;
using fluoroseq::hybrid_main;
using fluoroseq::print_bad_inputs;
using fluoroseq::print_invalid_classifier;
using fluoroseq::print_mpi_info;
}  // namespace

int main(int argc, char** argv) {
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif  // USE_MPI
    print_mpi_info();
    if (argc < 2) {
        print_bad_inputs();
        return 1;
    }
    char* mode = argv[1];
    int return_code;
    if (0 == strcmp(mode, "hmm")) {
        return_code = hmm_main(argc, argv);
    } else if (0 == strcmp(mode, "ann")) {
        return_code = ann_main(argc, argv);
    } else if (0 == strcmp(mode, "hybrid")) {
        return_code = hybrid_main(argc, argv);
    } else {
        print_invalid_classifier();
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif  // USE_MPI
    return return_code;
}
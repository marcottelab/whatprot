/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// For MPI version, define compiler macro USE_MPI when building.

// Standard C++ library headers:
#include <cstring>

// MPI header:
#ifdef USE_MPI
#include <mpi.h>
#endif  // USE_MPI

// Local project headers:
#include "main/classify/ann_main.h"
#include "main/classify/hmm_main.h"
#include "main/classify/hybrid_main.h"
#include "main/cmd_line_out.h"

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
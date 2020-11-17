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
#include "main/cmd_line_out.h"
#include "main/simulate/dt_main.h"
#include "main/simulate/rad_main.h"

namespace {
using fluoroseq::dt_main;
using fluoroseq::print_bad_inputs;
using fluoroseq::print_invalid_classifier;
using fluoroseq::print_mpi_info;
using fluoroseq::rad_main;
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
    if (0 == strcmp(mode, "rad")) {
        return_code = rad_main(argc, argv);
    } else if (0 == strcmp(mode, "dt")) {
        return_code = dt_main(argc, argv);
    } else {
        print_invalid_classifier();
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif  // USE_MPI
    return return_code;
}
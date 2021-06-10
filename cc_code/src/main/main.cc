/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Standard C++ library headers:
#include <cstring>

// Local project headers:
#include "main/classify/classify-main.h"
#include "main/cmd-line-out.h"
#include "main/fit/fit-main.h"
#include "main/simulate/simulate-main.h"

namespace {
using whatprot::classify_main;
using whatprot::fit_main;
using whatprot::print_bad_inputs;
using whatprot::print_invalid_command;
using whatprot::print_omp_info;
using whatprot::simulate_main;
}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_bad_inputs();
        return 1;
    }
    print_omp_info();
    char* mode = argv[1];
    int return_code;
    if (0 == strcmp(mode, "classify")) {
        return_code = classify_main(argc, argv);
    } else if (0 == strcmp(mode, "simulate")) {
        return_code = simulate_main(argc, argv);
    } else if (0 == strcmp(mode, "fit")) {
        return_code = fit_main(argc, argv);
    } else {
        return_code = 1;
        print_invalid_command();
    }
    return return_code;
}

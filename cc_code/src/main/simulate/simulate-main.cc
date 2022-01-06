/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Standard C++ library headers:
#include <cstring>
#include <iostream>

// Local project headers:
#include "main/cmd-line-out.h"
#include "main/simulate/dt-main.h"
#include "main/simulate/rad-main.h"

namespace whatprot {

namespace {
using std::cout;
}  // namespace

int simulate_main(int argc, char** argv) {
    if (argc < 3) {
        cout << "Usage: whatprot simulate [rad|dt]\n";
        return 1;
    }
    char* mode = argv[2];
    int return_code;
    if (0 == strcmp(mode, "rad")) {
        return_code = rad_main(argc, argv);
    } else if (0 == strcmp(mode, "dt")) {
        return_code = dt_main(argc, argv);
    } else {
        return_code = 1;
        print_invalid_classifier();
    }
    return return_code;
}

}  // namespace whatprot

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
#include "main/classify/hmm-main.h"
#include "main/classify/hybrid-main.h"
#include "main/classify/nn-main.h"
#include "main/cmd-line-out.h"

namespace whatprot {

int classify_main(int argc, char** argv) {
    if (argc < 3) {
        print_bad_inputs();
        return 1;
    }
    char* mode = argv[2];
    int return_code;
    if (0 == strcmp(mode, "hmm")) {
        return_code = hmm_main(argc, argv);
    } else if (0 == strcmp(mode, "hybrid")) {
        return_code = hybrid_main(argc, argv);
    } else if (0 == strcmp(mode, "nn")) {
        return_code = nn_main(argc, argv);
    } else {
        print_invalid_classifier();
    }
    return return_code;
}

}  // namespace whatprot

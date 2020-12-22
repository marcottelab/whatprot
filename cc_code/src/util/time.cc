/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "time.h"

#include <ctime>

namespace fluoroseq {

double wall_time() {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

double wall_tick() {
    return (double)1.0 / (double)CLOCKS_PER_SEC;
}

unsigned int time_based_seed() {
    return (unsigned int)clock();
}

}  // namespace fluoroseq
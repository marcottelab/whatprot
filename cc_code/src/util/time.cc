/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "time.h"

// Standard c++ headers
#include <ctime>

// OpenMP
#include <omp.h>

namespace whatprot {

double wall_time() {
    return omp_get_wtime();
}

double wall_tick() {
    return omp_get_wtick();
}

unsigned int time_based_seed() {
    return (unsigned int)clock();
}

}  // namespace whatprot
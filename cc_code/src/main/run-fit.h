/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_MAIN_RUN_FIT_H
#define WHATPROT_MAIN_RUN_FIT_H

// Standard C++ library headers:
#include <string>

namespace whatprot {

void run_fit(double stopping_threshold,
             double hmm_pruning_cutoff,
             std::string dye_seq_string,
             std::string radiometries_filename);

}  // namespace whatprot

#endif  // WHATPROT_MAIN_RUN_FIT_H

/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_MAIN_RUN_SIMULATE_DT_H
#define WHATPROT_MAIN_RUN_SIMULATE_DT_H

// Standard C++ library headers:
#include <string>

namespace whatprot {

void run_simulate_dt(unsigned int num_timesteps,
                     unsigned int num_to_generate,
                     std::string dye_seqs_filename,
                     std::string dye_tracks_filename);

}  // namespace whatprot

#endif  // WHATPROT_MAIN_RUN_SIMULATE_DT_H

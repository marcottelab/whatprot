/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_MAIN_RUN_SIMULATE_RAD_H
#define WHATPROT_MAIN_RUN_SIMULATE_RAD_H

#include <string>

namespace whatprot {

void run_simulate_rad(unsigned int num_timesteps,
                      unsigned int num_to_generate,
                      std::string dye_seqs_filename,
                      std::string radiometries_filename,
                      std::string ys_filename);

}  // namespace whatprot

#endif  // WHATPROT_MAIN_RUN_SIMULATE_RAD_H

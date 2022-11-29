/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_MAIN_RUN_CLASSIFY_NN_H
#define WHATPROT_MAIN_RUN_CLASSIFY_NN_H

#include <string>

namespace whatprot {

void run_classify_nn(std::string seq_params_filename,
                     int k,
                     double sig,
                     std::string dye_tracks_filename,
                     std::string radiometries_filename,
                     std::string predictions_filename);

}  // namespace whatprot

#endif  // WHATPROT_MAIN_RUN_CLASSIFY_NN_H
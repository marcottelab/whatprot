/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_MAIN_RUN_CLASSIFY_HMM_H
#define WHATPROT_MAIN_RUN_CLASSIFY_HMM_H

#include <string>

namespace whatprot {

void run_classify_hmm(std::string dye_seqs_filename,
                      std::string radiometries_filename,
                      std::string predictions_filename);

}  // namespace whatprot

#endif  // WHATPROT_MAIN_RUN_CLASSIFY_HMM_H
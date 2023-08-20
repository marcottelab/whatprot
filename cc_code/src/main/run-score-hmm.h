/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_MAIN_RUN_SCORE_HMM_H
#define WHATPROT_MAIN_RUN_SCORE_HMM_H

#include <string>

namespace whatprot {

void run_score_hmm(std::string seq_params_filename,
                   double hmm_pruning_cutoff,
                   std::string dye_seqs_filename,
                   std::string radiometries_filename,
                   std::string predictions_filename);

}  // namespace whatprot

#endif  // WHATPROT_MAIN_RUN_SCORE_HMM_H
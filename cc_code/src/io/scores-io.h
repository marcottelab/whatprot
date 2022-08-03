/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_IO_SCORES_IO_H
#define WHATPROT_IO_SCORES_IO_H

// Standard C++ library headers:
#include <string>
#include <vector>

namespace whatprot {

void write_scores(const std::string& filename,
                  const std::vector<std::vector<double>>& all_scores);

}  // namespace whatprot

#endif  // WHATPROT_IO_SCORES_IO_H
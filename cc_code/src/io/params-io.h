/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_IO_PARAMS_IO_H
#define WHATPROT_IO_PARAMS_IO_H

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

void write_params(const std::string& filename,
                  unsigned int num_channels,
                  const std::vector<SequencingModel>& models);

}  // namespace whatprot

#endif  // WHATPROT_IO_PARAMS_IO_H
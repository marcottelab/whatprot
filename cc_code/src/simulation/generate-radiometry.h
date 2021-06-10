/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_SIMULATION_GENERATE_RADIOMETRY_H
#define WHATPROT_SIMULATION_GENERATE_RADIOMETRY_H

// Standard C++ library headers:
#include <random>

// Local project headers:
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

void generate_radiometry(const SequencingModel& seq_model,
                         const DyeSeq& dye_seq,
                         unsigned int num_timesteps,
                         unsigned int num_channels,
                         std::default_random_engine* generator,
                         Radiometry* radiometry);

}  // namespace whatprot

#endif  // WHATPROT_SIMULATION_GENERATE_RADIOMETRY_H

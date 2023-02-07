/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_FITTERS_BOOTSTRAP_FIT_H
#define WHATPROT_FITTERS_BOOTSTRAP_FIT_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"

namespace whatprot {

void bootstrap_fit(unsigned int num_timesteps,
                   unsigned int num_channels,
                   double stopping_threshold,
                   const SequencingModel& seq_model,
                   const SequencingSettings& seq_settings,
                   const DyeSeq& dye_seq,
                   const std::vector<Radiometry>& radiometries,
                   unsigned int num_bootstrap_rounds,
                   double confidence_interval,
                   std::vector<SequencingModel>* seq_models);

}  // namespace whatprot

#endif  // WHATPROT_FITTERS_BOOTSTRAP_FIT_H

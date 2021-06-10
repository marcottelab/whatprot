/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "radiometry-precomputations.h"

// Local project headers:
#include "common/radiometry.h"
#include "hmm/step/peptide-emission.h"
#include "hmm/step/stuck-dye-emission.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

RadiometryPrecomputations::RadiometryPrecomputations(
        const Radiometry& radiometry,
        const SequencingModel& seq_model,
        int max_num_dyes)
        : peptide_emission(radiometry, max_num_dyes, seq_model) {
    for (unsigned int c = 0; c < radiometry.num_channels; c++) {
        stuck_dye_emissions.push_back(
                StuckDyeEmission(radiometry, c, seq_model));
    }
}

}  // namespace whatprot

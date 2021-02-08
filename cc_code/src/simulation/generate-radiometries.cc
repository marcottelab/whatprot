/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "generate-radiometries.h"

// Standard C++ library headers:
#include <random>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/error-model.h"
#include "common/radiometry.h"
#include "simulation/generate-radiometry.h"

namespace whatprot {

void generate_radiometries(
        const ErrorModel& error_model,
        const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        int num_timesteps,
        int num_channels,
        int radiometries_per_peptide,
        std::default_random_engine* generator,
        std::vector<SourcedData<Radiometry, SourceCount<int>>>* radiometries) {
    radiometries->reserve(dye_seqs.size() * radiometries_per_peptide);
    for (const SourcedData<DyeSeq, SourceCount<int>>& dye_seq : dye_seqs) {
        // We want to generate a certain number of radiometries per peptide,
        // not per dye_seq. Therefore we do this on repeat for each peptide
        // that produced this dye_seq.
        for (int i = 0; i < dye_seq.source.count; i++) {
            for (int j = 0; j < radiometries_per_peptide; j++) {
                radiometries->push_back(
                        SourcedData<Radiometry, SourceCount<int>>(
                                Radiometry(num_timesteps, num_channels),
                                dye_seq.source));
                generate_radiometry(error_model,
                                    dye_seq.value,
                                    num_timesteps,
                                    num_channels,
                                    generator,
                                    &radiometries->back().value);
                // Ignore any Radiometry with all 0s because it wouldn't be
                // detectable. Any Radiometry with all 0s at the 0th timestep
                // will have all 0s throughout.
                bool nontrivial = false;
                for (int c = 0; c < num_channels; c++) {
                    if (radiometries->back().value(0, c) != 0.0) {
                        nontrivial = true;
                    }
                }
                if (!nontrivial) {
                    radiometries->pop_back();
                }
            }
        }
    }
}

}  // namespace whatprot

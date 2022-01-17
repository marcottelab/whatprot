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
#include "common/radiometry.h"
#include "parameterization/model/sequencing-model.h"
#include "simulation/generate-radiometry.h"

namespace whatprot {

namespace {
using std::default_random_engine;
using std::discrete_distribution;
using std::vector;
}  // namespace

void generate_radiometries(
        const SequencingModel& seq_model,
        const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        unsigned int num_timesteps,
        unsigned int num_channels,
        unsigned int num_to_generate,
        std::default_random_engine* generator,
        std::vector<SourcedData<Radiometry, SourceCount<int>>>* radiometries) {
    // We want the dye tracks generated based on a uniform distribution of
    // peptides, not radiometries. We therefore need a discrete_distribution,
    // because it is weighted.
    vector<double> index_to_weight;
    for (unsigned int i = 0; i < dye_seqs.size(); i++) {
        index_to_weight[i] = (double)dye_seqs[i].source.count;
    }
    discrete_distribution<unsigned int> random_dye_seq_idx(
            index_to_weight.begin(), index_to_weight.end());
    for (unsigned int i = 0; i < num_to_generate; i++) {
        unsigned int dye_seq_idx = random_dye_seq_idx(*generator);
        radiometries->push_back(SourcedData<Radiometry, SourceCount<int>>(
                Radiometry(num_timesteps, num_channels),
                dye_seqs[dye_seq_idx].source));
        generate_radiometry(seq_model,
                            dye_seqs[dye_seq_idx].value,
                            num_timesteps,
                            num_channels,
                            generator,
                            &radiometries->back().value);
        // Ignore any Radiometry with all 0s because it wouldn't be detectable.
        // Any Radiometry with all 0s at the 0th timestep will have all 0s
        // throughout.
        bool trivial = true;
        for (unsigned int c = 0; c < num_channels; c++) {
            if (radiometries->back().value(0, c) != 0.0) {
                trivial = false;
            }
        }
        if (trivial) {
            radiometries->pop_back();
        }
    }
}

}  // namespace whatprot
